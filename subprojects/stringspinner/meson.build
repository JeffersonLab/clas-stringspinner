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

# dependency pythia8
pythia_link_args = []
pythia_config_prog = find_program('pythia8-config', required: false)
if not pythia_config_prog.found()
  error('Pythia8 is not found')
endif
foreach arg : run_command(pythia_config_prog, '--libs', check: true).stdout().split()
  if not arg.contains('-rpath')
    pythia_link_args += arg
  endif
endforeach

# build stringspinner libraries
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

# build stringspinner example executable
stringspinner_inc = include_directories('stringspinner')
dis_exe = executable(
  'stringspinner-example',
  sources: ['stringspinner/dis.cc'],
  include_directories: stringspinner_inc,
  link_with: [ def_lib, mc3P0_lib ],
  link_args: pythia_link_args,
  install: true,
  cpp_args: ['-Wno-sign-compare', '-Wno-unused-variable'],
)

# declare subproject
stringspinner_dep = declare_dependency(
  include_directories: stringspinner_inc,
  link_with: [ def_lib, mc3P0_lib ],
)
