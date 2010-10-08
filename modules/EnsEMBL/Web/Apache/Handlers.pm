package EnsEMBL::Web::Apache::Handlers;

use strict;

use Apache2::Const qw(:common :http :methods);
use Apache2::SizeLimit;
use Apache2::Connection;
use Apache2::URI;
use APR::URI;
use CGI::Cookie;
use Config;
use Fcntl ':flock';
use Sys::Hostname;
use Time::HiRes qw(time);
use URI::Escape qw(uri_escape);

use SiteDefs qw(:APACHE);

use Bio::EnsEMBL::Registry;

use EnsEMBL::Web::Cache;
use EnsEMBL::Web::Cookie;
use EnsEMBL::Web::Registry;
use EnsEMBL::Web::RegObj;
use EnsEMBL::Web::SpeciesDefs;

use EnsEMBL::Web::Apache::DasHandler;
use EnsEMBL::Web::Apache::SSI;
use EnsEMBL::Web::Apache::SpeciesHandler;

our $species_defs = new EnsEMBL::Web::SpeciesDefs;
our $MEMD         = new EnsEMBL::Web::Cache;

our $BLAST_LAST_RUN;
our $LOAD_COMMAND;

BEGIN {
  $LOAD_COMMAND = $Config{'osname'} eq 'dec_osf' ? \&_load_command_alpha :
                  $Config{'osname'} eq 'linux'   ? \&_load_command_linux :
                                                   \&_load_command_null;
};

#======================================================================#
# Perl apache handlers in order they get executed                      #
#======================================================================#

# Child Init Handler
# Sets up the web registry object - and initializes the timer
sub childInitHandler {
  my $r = shift;
  
  my @X             = localtime;
  my $temp_hostname = hostname;
  my $temp_proc_id  = '' . reverse $$;
  my $temp_seed     = ($temp_proc_id + $temp_proc_id << 15) & 0xffffffff;
  
  while ($temp_hostname =~ s/(.{1,4})//) {
    $temp_seed = $temp_seed ^ unpack("%32L*", $1);
  }
  
  srand(time ^ $temp_seed);
  
  # Create the Registry
  $ENSEMBL_WEB_REGISTRY = new EnsEMBL::Web::Registry;
  $ENSEMBL_WEB_REGISTRY->timer->set_process_child_count(0);
  $ENSEMBL_WEB_REGISTRY->timer->set_process_start_time(time);
  
  warn sprintf "Child initialised: %7d %04d-%02d-%02d %02d:%02d:%02d\n", $$, $X[5]+1900, $X[4]+1, $X[3], $X[2], $X[1], $X[0] if $ENSEMBL_DEBUG_FLAGS & $SiteDefs::ENSEMBL_DEBUG_HANDLER_ERRORS;
}

sub redirect_to_nearest_mirror {
  my $r = shift;
  
  if ($species_defs->ENSEMBL_MIRRORS) {
    my $unparsed_uri = $r->unparsed_uri;

    ## Check url
    if ($unparsed_uri =~ /redirect=mirror/) {
      my ($referer) = $unparsed_uri =~ /source=([\w\.-]+)/;
      
      ## Display the redirect message (but only if user comes from other mirror)
      if ($referer && $referer ne $species_defs->ENSEMBL_SERVERNAME && grep { $referer =~ /$_/ } values %{$species_defs->ENSEMBL_MIRRORS}) {
        my $back = 'http://' . $referer . $unparsed_uri;
        $back =~ s/;?source=$referer//;

        my $user_message = qq{You've been redirected to your nearest mirror - } . $species_defs->ENSEMBL_SERVERNAME . "\n";
        $user_message   .= qq{<ul><li>Take me back to <a href="$back">$referer</a></li></ul>};
        
        my $cookie = new CGI::Cookie(
          -name    => 'user_message',
          -value   => uri_escape($user_message),
          -expires => '+1m',     
        );
        
        ## Redirecting to same page, but without redirect params in url
        $r->err_headers_out->add('Set-Cookie' => $cookie);
        
        $unparsed_uri =~ s/;?source=$referer//;
        $unparsed_uri =~ s/;?redirect=mirror//;
        $unparsed_uri =~ s/\?$//;
        
        $r->headers_out->set(Location => $species_defs->ENSEMBL_BASE_URL . $unparsed_uri);
        
        return REDIRECT;
      }
      
      my $cookie = new CGI::Cookie(
        -name    => 'redirect',
        -value   => 'mirror',
        -expires => 'Thu, 31-Dec-2037 22:22:22 GMT', ## End of time :)     
      );
      
      $r->err_headers_out->add('Set-Cookie' => $cookie);

      return DECLINED;
    }

    ## Check "don't redirect me" cookie
    my %cookies = CGI::Cookie->parse($r->headers_in->{'Cookie'});
    
    return DECLINED if $cookies{'redirect'} && $cookies{'redirect'}->value eq 'mirror';
    
    if ($species_defs->ENSEMBL_MIRRORS && keys %{$species_defs->ENSEMBL_MIRRORS}) {
      my $geo;
      
      eval '
        use Geo::IP;
        $geo = new Geo::IP(GEOIP_MEMORY_CACHE | GEOIP_CHECK_CACHE);
      ';
      
      warn $@ if $@;
    
      ## Ok, so which country you from
      if ($geo) {
        my $ip       = $r->headers_in->{'X-Forwarded-For'} || $r->connection->remote_ip;
        my $country  = $geo->country_code_by_addr($ip);
        my $location = $species_defs->ENSEMBL_MIRRORS->{$country} || $species_defs->ENSEMBL_MIRRORS->{'MAIN'};
        
        if ($location) {
          return DECLINED if $location eq $species_defs->ENSEMBL_SERVERNAME;
    
          ## Deleting cookie for current site
          my $cookie = new CGI::Cookie(
            -name    => 'redirect',
            -value   => '',
            -expires => '-1h',         
          );
    
          $unparsed_uri .= $unparsed_uri =~ /\?/ ? ';redirect=mirror' : '?redirect=mirror';
          $unparsed_uri .= ';source=' . $species_defs->ENSEMBL_SERVERNAME;
          
          $r->err_headers_out->add('Set-Cookie' => $cookie);
          $r->headers_out->set(Location => "http://$location$unparsed_uri");
          
          return REDIRECT;       
        }
      }
    }
  }

  return DECLINED;
}

sub postReadRequestHandler {
  my $r = shift; # Get the connection handler

  # Nullify tags
  $ENV{'CACHE_TAGS'} = {};
  
  # Manipulate the Registry
  $ENSEMBL_WEB_REGISTRY->timer->new_child;
  $ENSEMBL_WEB_REGISTRY->timer->clear_times;
  $ENSEMBL_WEB_REGISTRY->timer_push('Handling script', undef, 'Apache');
  
  my $user_cookie = new EnsEMBL::Web::Cookie({
    host    => $ENSEMBL_COOKIEHOST,
    name    => $ENSEMBL_USER_COOKIE,
    value   => '',
    env     => 'ENSEMBL_USER_ID',
    hash    => {
      offset  => $ENSEMBL_ENCRYPT_0,
      key1    => $ENSEMBL_ENCRYPT_1,
      key2    => $ENSEMBL_ENCRYPT_2,
      key3    => $ENSEMBL_ENCRYPT_3,
      expiry  => $ENSEMBL_ENCRYPT_EXPIRY,
      refresh => $ENSEMBL_ENCRYPT_REFRESH
    }
  });
  
  # Initialize the user (and possibly group) objects
  # Unlikely to go to db - just store the IDs
  $ENSEMBL_WEB_REGISTRY->initialize_user($user_cookie, $r);
  
  ## Ajax cookie
  my %cookies = CGI::Cookie->parse($r->headers_in->{'Cookie'});
  
  $r->subprocess_env->{'ENSEMBL_AJAX_VALUE'}  = $cookies{'ENSEMBL_AJAX'}  && $cookies{'ENSEMBL_AJAX'}->value  ? $cookies{'ENSEMBL_AJAX'}->value  : 'none';
  $r->subprocess_env->{'ENSEMBL_IMAGE_WIDTH'} = $cookies{'ENSEMBL_WIDTH'} && $cookies{'ENSEMBL_WIDTH'}->value ? $cookies{'ENSEMBL_WIDTH'}->value : ($ENSEMBL_IMAGE_WIDTH || 800);
  
  $ENSEMBL_WEB_REGISTRY->timer_push('Post read request handler completed', undef, 'Apache');
  
  # Ensembl DEBUG cookie
  $r->headers_out->add('X-MACHINE' => $ENSEMBL_SERVER) if $cookies{'ENSEMBL_DEBUG'};
  
  return;
}

sub cleanURI {
  my $r = shift;
  
  # Void call to populate ENV
  $r->subprocess_env;
  
  # Clean out the uri
  my $uri = $ENV{'REQUEST_URI'};
  
  if ($uri =~ s/time=\d+\.\d+[;&]?//g + $uri =~ s!([^:])/{2,}!$1/!g) {
    $r->parse_uri($uri);
    $r->subprocess_env->{'REQUEST_URI'} = $uri;
  }

  # Clean out the referrer
  my $referer = $ENV{'HTTP_REFERER'};
  
  if ($referer =~ s/time=\d+\.\d+[;&]?//g + $referer =~ s!([^:])/{2,}!$1/!g) {
    $r->subprocess_env->{'HTTP_REFERER'} = $referer;
  }
  
  return DECLINED;
}

sub handler {
  my $r = shift; # Get the connection handler
  
  $ENSEMBL_WEB_REGISTRY->timer->set_name('REQUEST ' . $r->uri);
  
  my $u           = $r->parsed_uri;
  my $file        = $u->path;
  my $querystring = $u->query;
  
  my $session_cookie = new EnsEMBL::Web::Cookie({
    host    => $ENSEMBL_COOKIEHOST,
    name    => $ENSEMBL_SESSION_COOKIE,
    value   => '',
    env     => 'ENSEMBL_SESSION_ID',
    hash    => {
      offset  => $ENSEMBL_ENCRYPT_0,
      key1    => $ENSEMBL_ENCRYPT_1,
      key2    => $ENSEMBL_ENCRYPT_2,
      key3    => $ENSEMBL_ENCRYPT_3,
      expiry  => $ENSEMBL_ENCRYPT_EXPIRY,
      refresh => $ENSEMBL_ENCRYPT_REFRESH
    }
  });

  my @raw_path = split m|/|, $file;
  shift @raw_path; # Always empty

  my %species_lookup = map { $_ => 1 } $species_defs->valid_species;
  
  my %species_map = (
    common => 'common',
    multi  => 'Multi',
    perl   => $ENSEMBL_PRIMARY_SPECIES,
    map { lc($_) => $ENSEMBL_SPECIES_ALIASES->{$_} } keys %$ENSEMBL_SPECIES_ALIASES
  );

  $species_map{lc $_} = $_ for values %species_map; # Self-mapping
  
  ## Identify the species element, if any
  my ($species, @path_segments);
 
  ## Check for stable id URL (/id/ENSG000000nnnnnn)
  if ($raw_path[0] && $raw_path[0] =~ /^id$/i && $raw_path[1]) {
    my $stable_id = $raw_path[1];
    my ($object_type, $db_type, $uri);
    ($species, $object_type, $db_type) = Bio::EnsEMBL::Registry->get_species_and_object_type($stable_id);
    
    if ($species && $object_type) {
      $uri = "/$species/";
      
      if ($object_type eq 'Gene') {
        $uri .= "Gene/Summary?g=$stable_id";
      } elsif ($object_type eq 'Transcript') {
        $uri .= "Transcript/Summary?t=$stable_id";
      } elsif ($object_type eq 'Translation') {
        $uri .= "Transcript/ProteinSummary?t=$stable_id";
      } else {
        $uri .= "psychic?q=$stable_id";
      }
      
      $r->uri($uri);
      $r->headers_out->add('Location' => $r->uri);
      $r->child_terminate;
      
      $ENSEMBL_WEB_REGISTRY->timer_push('Handler "REDIRECT"', undef, 'Apache');
    
      return HTTP_MOVED_PERMANENTLY;
    }
  }
  
  foreach (@raw_path) {
    if ($species_lookup{$species_map{lc $_}} && !$species) {
      $species = $_;
    } else {
      push @path_segments, $_;
    }
  }
  
  if (!$species) {
    if (grep /$raw_path[0]/, qw(Multi das common default)) {
      $species = $raw_path[0];
      shift @path_segments;
    } elsif ($path_segments[0] eq 'Gene' && $querystring) {
      my %param = split ';|=', $querystring;
      if (my $gene_stable_id = $param{'g'}) {
        my ($id_species) = Bio::EnsEMBL::Registry->get_species_and_object_type($gene_stable_id);
        
        $species = $id_species if $id_species;
      }  
    }
  }
 
  @path_segments = @raw_path unless $species;
  
  # Some memcached tags (mainly for statistics)
  my $prefix = '';
  my @tags = map { $prefix = join '/', $prefix, $_; $prefix; } @path_segments;
  @tags = map {( "/$species$_", $_ )} @tags;
  push @tags, "/$species";
  
  $ENV{'CACHE_TAGS'}{$_} = 1 for @tags; # /memcached tags
    
  my $Tspecies  = $species;
  my $script    = undef;
  my $path_info = undef;
  my $species_name = $species_map{lc $species};
  my $return;
  
  if (!$species && $raw_path[-1] !~ /\./) {
    $species      = 'common';
    $species_name = 'common';
    $file         = "/common$file";
    $file         =~ s|/$||;
  }
  
  if ($raw_path[0] eq 'das') {
    my ($das_species) = split /\./, $path_segments[0];
    
    $return = EnsEMBL::Web::Apache::DasHandler::handler_das($r, $session_cookie, $species_map{lc $das_species}, \@path_segments, $querystring);
    
    $ENSEMBL_WEB_REGISTRY->timer_push('Handler for DAS scripts finished', undef, 'Apache');
  } elsif ($species && $species_name) { # species script
    $return = EnsEMBL::Web::Apache::SpeciesHandler::handler_species($r, $session_cookie, $species_name, \@path_segments, $querystring, $file, $species_name eq $species);
    
    $ENSEMBL_WEB_REGISTRY->timer_push('Handler for species scripts finished', undef, 'Apache');
    
    shift @path_segments;
    shift @path_segments;
  }
  
  if (defined $return) {
    if ($return == OK) {
      push_script_line($r) if $ENSEMBL_DEBUG_FLAGS & $SiteDefs::ENSEMBL_DEBUG_HANDLER_ERRORS;
      
      $r->push_handlers(PerlCleanupHandler => \&cleanupHandler_script);
      $r->push_handlers(PerlCleanupHandler => \&Apache2::SizeLimit::handler);
    }
    
    return $return;
  }
  
  $species = $Tspecies;
  $script = join '/', @path_segments;

  # Permanent redirect for old species home pages:
  # e.g. /Homo_sapiens or Homo_sapiens/index.html -> /Homo_sapiens/Info/Index
  if ($species && $species_name && (!$script || $script eq 'index.html')) {
    $r->uri($species_name eq 'common' ? 'index.html' : "/$species_name/Info/Index");
    $r->headers_out->add('Location' => $r->uri);
    $r->child_terminate;
    $ENSEMBL_WEB_REGISTRY->timer_push('Handler "REDIRECT"', undef, 'Apache');
    
    return HTTP_MOVED_PERMANENTLY;
  }

  return DECLINED if $species eq 'biomart' && $script =~ /^mart(service|results|view)/;

  my $path = join '/', $species || (), $script || (), $path_info || ();
  
  $r->uri("/$path");
  
  my $filename = $MEMD ? $MEMD->get("::STATIC::$path") : '';
  
  # Search the htdocs dirs for a file to return
  # Exclude static files (and no, html is not a static file in ensembl)
  if ($path !~ /\.(\w{2,3})$/) {
    if (!$filename) {
      foreach my $dir (grep { -d $_ && -r $_ } @ENSEMBL_HTDOCS_DIRS) {
        my $f = "$dir/$path";
        
        if (-d $f || -r $f) {
          $filename = -d $f ? '! ' . $f : $f;
          $MEMD->set("::STATIC::$path", $filename, undef, 'STATIC') if $MEMD;
          
          last;
        }
      }
    }
  }
  
  if ($filename =~ /^! (.*)$/) {
    $r->uri($r->uri . ($r->uri      =~ /\/$/ ? '' : '/') . 'index.html');
    $r->filename($1 . ($r->filename =~ /\/$/ ? '' : '/') . 'index.html');
    $r->headers_out->add('Location' => $r->uri);
    $r->child_terminate;
    $ENSEMBL_WEB_REGISTRY->timer_push('Handler "REDIRECT"', undef, 'Apache');
    
    return HTTP_TEMPORARY_REDIRECT;
  } elsif ($filename) {
    $r->filename($filename);
    $r->content_type('text/html');
    $ENSEMBL_WEB_REGISTRY->timer_push('Handler "OK"', undef, 'Apache');
    
    EnsEMBL::Web::Apache::SSI::handler($r);
    
    return OK;
  }
  
  # Give up
  $ENSEMBL_WEB_REGISTRY->timer_push('Handler "DECLINED"', undef, 'Apache');
  
  return DECLINED;
}

sub logHandler {
  my $r = shift;
  my $T = time;
  
  $r->subprocess_env->{'ENSEMBL_CHILD_COUNT'}  = $ENSEMBL_WEB_REGISTRY->timer->get_process_child_count;
  $r->subprocess_env->{'ENSEMBL_SCRIPT_START'} = sprintf '%0.6f', $T;
  $r->subprocess_env->{'ENSEMBL_SCRIPT_END'}   = sprintf '%0.6f', $ENSEMBL_WEB_REGISTRY->timer->get_script_start_time;
  $r->subprocess_env->{'ENSEMBL_SCRIPT_TIME'}  = sprintf '%0.6f', $T - $ENSEMBL_WEB_REGISTRY->timer->get_script_start_time;
  
  return DECLINED;
}

sub cleanupHandler {
  my $r = shift;  # Get the connection handler
  
  return if $r->subprocess_env->{'ENSEMBL_ENDTIME'};
  
  my $end_time   = time;
  my $start_time = $ENSEMBL_WEB_REGISTRY->timer->get_script_start_time;
  my $length     = $end_time - $start_time;
  
  if ($length >= $ENSEMBL_LONGPROCESS_MINTIME) {
    my $u     = $r->parsed_uri;
    my $file  = $u->path;
    my $query = $u->query . $r->subprocess_env->{'ENSEMBL_REQUEST'};
    my $size  = &$Apache2::SizeLimit::HOW_BIG_IS_IT();
    
    $r->subprocess_env->{'ENSEMBL_ENDTIME'} = $end_time;
    
    if ($ENSEMBL_DEBUG_FLAGS & $SiteDefs::ENSEMBL_DEBUG_HANDLER_ERRORS) {
      my @X = localtime($start_time);
      
      warn sprintf(
        "LONG PROCESS: %12s DT:  %04d-%02d-%02d %02d:%02d:%02d Time: %10s Size: %10s\nLONG PROCESS: %12s REQ: %s\nLONG PROCESS: %12s IP:  %s  UA: %s\n", 
        $$, $X[5]+1900, $X[4]+1, $X[3], $X[2], $X[1], $X[0], $length, $size, 
        $$, "$file?$query", 
        $$, $r->subprocess_env->{'HTTP_X_FORWARDED_FOR'}, $r->headers_in->{'User-Agent'}
      );
    }
  }

  # Now we check if the die file has been touched...
  my $die_file = $ENSEMBL_SERVERROOT . '/logs/ensembl.die';
  
  if (-e $die_file) {
    my @temp = stat $die_file;
    my $file_mod_time = $temp[9];
    if ($file_mod_time >= $ENSEMBL_WEB_REGISTRY->timer->get_process_start_time) {
      warn sprintf "KILLING CHILD %10s\n", $$;
      
      if ($Apache2::SizeLimit::WIN32) {
        CORE::exit(-2);
      } else {
        $r->child_terminate;
      }
    }
    
    return DECLINED;
  }
}

sub cleanupHandler_script {
  my $r = shift;
  
  $ENSEMBL_WEB_REGISTRY->timer_push('Cleaned up', undef, 'Cleanup');
  
  warn $ENSEMBL_WEB_REGISTRY->timer->render if $ENSEMBL_DEBUG_FLAGS & $SiteDefs::ENSEMBL_DEBUG_PERL_PROFILER;
  
  push_script_line($r, 'ENDSCR', sprintf '%10.3f', time - $r->subprocess_env->{'LOG_TIME'}) if $ENSEMBL_DEBUG_FLAGS & $SiteDefs::ENSEMBL_DEBUG_HANDLER_ERRORS;
  
  cleanupHandler_blast($r) if $ENSEMBL_BLASTSCRIPT;
}

sub cleanupHandler_blast {
  my $r = shift;
  
  my $directory = $ENSEMBL_TMP_DIR_BLAST . '/pending';
  my $FLAG  = 0;
  my $count = 0;
  my $ticket;
  my $_process_blast_called_at = time;

  $ticket = $ENV{'ticket'};
  
  # Lets work out when to run this!
  my $run_blast;
  my $loads = _get_loads();
  my $seconds_since_last_run = (time - $BLAST_LAST_RUN);

  if ($ticket) {
    if (_run_blast_ticket($loads, $seconds_since_last_run)) {
      $FLAG = 1;
      $BLAST_LAST_RUN = time;
    }
  } else {
    # Current run blasts
    if (_run_blast_no_ticket($loads, $seconds_since_last_run)) {
      $BLAST_LAST_RUN = time;
      $FLAG = 1;
    }
  }
  
  while ($FLAG) {
    $count++;
    $FLAG = 0;
    
    if (opendir(DH, $directory)) {
      while (my $FN = readdir(DH)) {
        my $file = "$directory/$FN";
        
        next unless -f $file; # File
        next if -z $file;     # Contains something
        
        my @STAT = stat $file;
        
        next if $STAT[8]+5 > time; # Was last modified more than 5 seconds ago
        next if $ticket && $file !~ /$ticket/;
        
        # We have a ticket
        open  FH, $file;
        
        flock FH, LOCK_EX;
        my $blast_file = <FH>;
        chomp $blast_file;
        
        $blast_file = $1 if $blast_file =~ /^([\/\w\.-]+)/;
        
        (my $FILE2 = $file) =~ s/pending/parsing/;
        
        rename $file, $FILE2;
        
        (my $FILE3 = $file) =~ s/pending/sent/;
        
        unlink $FILE3;
        
        flock FH, LOCK_UN;
        
        my $COMMAND = "$ENSEMBL_BLASTSCRIPT $blast_file $FILE2";
        
        warn "BLAST: $COMMAND";
        
        `$COMMAND`; # Now we parse the blast file
        
        if ($ticket && ($_process_blast_called_at + 30 > time)) {
          $loads = _get_loads();
          $FLAG = 1 if $count < 15;
        }
        
        last;
      }
      
      closedir(DH);
    }
  }
}

sub childExitHandler {
  my $r = shift;
  
  $ENSEMBL_WEB_REGISTRY->tidy_up if $ENSEMBL_WEB_REGISTRY; # Disconnect from the DB
  
  if ($ENSEMBL_DEBUG_FLAGS & $SiteDefs::ENSEMBL_DEBUG_HANDLER_ERRORS) {
    warn sprintf "Child %9d: - reaped at      %30s;  Time: %11.6f;  Req:  %4d;  Size: %8dK\n",
      $$, '' . gmtime, time-$ENSEMBL_WEB_REGISTRY->timer->get_process_start_time,
      $ENSEMBL_WEB_REGISTRY->timer->get_process_child_count,
      &$Apache2::SizeLimit::HOW_BIG_IS_IT()
  }
}

sub push_script_line {
  my $r      = shift;
  my $prefix = shift || 'SCRIPT';
  my $extra  = shift;
  my @X      = localtime;
  
  warn sprintf(
    "%s: %s%9d %04d-%02d-%02d %02d:%02d:%02d %s %s\n",
    $prefix, hostname, $$,
    $X[5] + 1900, $X[4] + 1, $X[3], $X[2], $X[1], $X[0],
    $r->subprocess_env->{'REQUEST_URI'}, $extra
  );
  
  $r->subprocess_env->{'LOG_TIME'} = time;
}

#======================================================================#
# BLAST Support functionality - TODO: update before implementing!      #
#======================================================================#

sub _run_blast_no_ticket {
  my ($loads, $seconds_since_last_run) = @_;
  return $loads->{'blast'} < 3 && rand $loads->{'httpd'} < 10 && rand $seconds_since_last_run > 1;
}

sub _run_blast_ticket {
  my ($loads, $seconds_since_last_run) = @_;
  return $loads->{'blast'} < 8;
}

sub _get_loads {
  return {
    blast => &$LOAD_COMMAND('blast'),
    httpd => &$LOAD_COMMAND('httpd')
  };
}

sub  _load_command_null {
  return 1;
}

sub _load_command_alpha {
  my $command = shift;
  my $VAL = `ps -A | grep $command | wc -l`;
  
  return $VAL - 1;
}

sub _load_command_linux {
  my $command = shift;
  my $VAL = `ps --no-heading -C $command  | wc -l`;
  
  return $VAL + 0;
}

1;
