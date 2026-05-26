project(
  'clas-stringspinner',
  'cpp', 'fortran'
  meson_version: '>=1.3',
  default_options: {
    'cpp_std':    'c++20',
    'buildtype':  'release',
    'libdir':     'lib',
  },
  version: '0.1.6',
)

# dependencies
fmt_dep    = dependency('fmt', method: 'pkg-config', version:                                       '>= 9.1.0')
pythia_dep = dependency('pythia', fallback: ['stringspinner', 'stringspinner_pythia_dep'], version: '>= 8.312')
hipo_dep   = dependency('hipo4', method: 'pkg-config', required: false, version:                    '>= 4.2.0')

# preprocessor macros
add_project_arguments('-D CLAS_STRINGSPINNER_VERSION="' + meson.project_version() + '"', language: ['cpp'])
if hipo_dep.found()
  add_project_arguments('-D CLAS_STRINGSPINNER_HIPO', language: ['cpp'])
endif

# stringspinner: Comput.Phys.Commun. 292 (2023) 108886; arXiv: https://arxiv.org/abs/2305.05058
install_headers(
  [
    'stringspinner/PrimordialKT.h',
    'stringspinner/StringSpinner.h',
    'stringspinner/Transversity.h',
    'stringspinner/VectorMesonDecays.h',
  ],
  subdir: 'stringspinner'
)
stringspinner_def_lib = static_library(
  'StringSpinner_def',
  'stringspinner/definitions.f90',
  install: true,
  fortran_args: ['-Wno-tabs'],
)
stringspinner_mc3P0_lib = static_library(
  'StringSpinner_mc3P0',
  'stringspinner/mc3P0.f90',
  install: true,
  link_with: stringspinner_def_lib,
  fortran_args: ['-frecord-marker=8', '-fbounds-check', '-Wno-tabs', '-Wno-unused-variable', '-Wno-conversion', '-Wno-unused-dummy-argument', '-Wno-aliasing'],
)
stringspinner_dep = declare_dependency(
  include_directories: include_directories('stringspinner', is_system: true), # use `is_system` to suppress warnings in consumer builds
  link_with: [stringspinner_def_lib, stringspinner_mc3P0_lib],
  dependencies: [pythia_dep],
)
if get_option('install_stringspinner_example')
  dis_exe = executable(
    'stringspinner-example',
    sources: ['stringspinner/dis.cc'],
    dependencies: [ stringspinner_dep, pythia_dep ],
    install: true,
    cpp_args: ['-Wno-sign-compare', '-Wno-unused-variable'],
  )
endif

# main executable
main_exe = executable(
  meson.project_name(),
  'src' / meson.project_name() + '.cpp',
  'src' / 'config' / 'src' / 'common.cc',
  dependencies: [ fmt_dep, stringspinner_dep, pythia_dep, hipo_dep ],
  install: true,
)

# test executable
test(
  meson.project_name(),
  main_exe,
  args: [
    '--num-events', '10000',
    '--pol-type', 'LU',
    '--beam-spin', 'p',
    '--cut-inclusive', '11,211,-211',
    '--seed', '828203788',
  ]
)

# install scripts
install_data(
  'src' / 'clas-stringspinner-slurm.rb',
  rename: 'clas-stringspinner-slurm',
  install_dir: get_option('bindir'),
  install_mode: 'rwxr-xr-x',
)
