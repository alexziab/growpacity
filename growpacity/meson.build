growpacity_extension = py3.extension_module(
  'growpacity_c',                # name of the module
  ['growpacity.c', 'growpacity_python.c'], # C-files
  include_directories: inc_np, # what directories to include
  install: true,
)

executable(
  'growpacity_cli',           # command line name
  'growpacity.c',             # needs to have the main function inside
  include_directories: inc_np,
  link_args: ['-lm'],
  install: true,
)

py3.install_sources(['__init__.py', 'growpacity.py'], subdir: 'growpacity')

subdir('data')
