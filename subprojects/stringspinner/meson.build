project(
  'stringspinner',
  'cpp', 'fortran',
  meson_version: '>=1.2',
  default_options: {
    'cpp_std':   'c++11',
    'buildtype': 'release',
    'libdir':    'lib',
  },
  version: '0.0.0',
)

pythia_link_args = []
foreach arg : run_command('pythia8-config', '--libs', check: true).stdout().split()
  if not arg.contains('-rpath')
    pythia_link_args += arg
  endif
endforeach

def_lib = static_library(
  'StringSpinner_def',
  'stringspinner/definitions.f90',
  install: true,
  fortran_args: ['-Wno-tabs'],
)

mc3P0_lib = static_library(
  'StringSpinner_mc3P0',
  'stringspinner/mc3P0.f90',
  install: true,
  link_with: def_lib,
  fortran_args: ['-frecord-marker=8', '-fbounds-check', '-Wno-tabs', '-Wno-unused-variable', '-Wno-conversion', '-Wno-unused-dummy-argument'],
)

dis_exe = executable(
  'dis',
  sources: ['stringspinner/dis.cc'],
  include_directories: include_directories('stringspinner'),
  link_with: [ def_lib, mc3P0_lib ],
  link_args: pythia_link_args,
  install: true,
  cpp_args: ['-Wno-sign-compare', '-Wno-unused-variable'],
)
