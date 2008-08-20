package EnsEMBL::Web::Component::Account::Login;

### Module to create user login form 

use strict;
use warnings;
no warnings "uninitialized";
use base qw(EnsEMBL::Web::Component::Account);
use EnsEMBL::Web::Form;

sub _init {
  my $self = shift;
  $self->cacheable( 0 );
  $self->ajaxable(  0 );
}

sub caption {
  my $self = shift;
  return 'Login';
}

sub content {
  my $self = shift;

  my $referer = CGI::escape($self->object->param('_referer'));
  my $form = EnsEMBL::Web::Form->new( 'login', "/Account/SetCookie", 'post' );
  my $reg_url = $self->url("/Account/Register?_referer=$referer");
  my $pwd_url = $self->url("/Account/LostPassword?_referer=$referer");

  $form->add_element('type'  => 'Email',    'name'  => 'email', 'label' => 'Email', 'required' => 'yes');
  $form->add_element('type'  => 'Password', 'name'  => 'password', 'label' => 'Password', 'required' => 'yes');
  $form->add_element('type'  => 'Hidden',   'name'  => 'url', 'value' => $self->object->param('url'));
  $form->add_element('type'  => 'Hidden',   'name'  => '_referer', 'value' => $referer);
  $form->add_element('type'  => 'Submit',   'name'  => 'submit', 'value' => 'Log in', 'class'=>'cp-refresh');
  $form->add_element('type'  => 'Information',
                     'value' => qq(<p><a href="$reg_url" class="cp-internal">Register</a>
                                  | <a href="$pwd_url" class="cp-internal">Lost password</a></p>));

  return $form->render;
}

1;
