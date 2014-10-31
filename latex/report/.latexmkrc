add_cus_dep('glo', 'gls', 0, 'makeglo2gls');
add_cus_dep('acn', 'acr', 0, 'makeglo2gls');
sub makeglo2gls {
  my ($base_name, $path) = fileparse( $_[0] );
  pushd($path);
  system("makeglossaries $base_name");
  popd;
}

add_cus_dep('nlo', 'nls', 0, 'makenlo2nls');
sub makenlo2nls {
  my ($base_name, $path) = fileparse( $_[0] );
  pushd($path);
  system( "makeindex -s nomencl.ist -o \"$base_name.nls\" \"$base_name.nlo\"" );
  popd;
}
