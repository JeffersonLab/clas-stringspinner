project(
  'clas-stringspinner',
  'cpp',
  meson_version: '>=1.3',
  default_options: {
    'cpp_std':    'c++20',
    'buildtype':  'release',
    'libdir':     'lib',
    'wrap_mode':  'forcefallback',
  },
  version: '0.1.6',
)

# subprojects
stringspinner_proj = subproject('stringspinner')

# dependencies
fmt_dep           = dependency('fmt', method: 'pkg-config', version: '>= 9.1.0')
stringspinner_dep = dependency('stringspinner', fallback: ['stringspinner', 'stringspinner_dep'])
pythia_dep        = dependency('pythia', fallback: ['stringspinner', 'stringspinner_pythia_dep'], version: '8.312')
hipo_dep          = dependency('hipo4', method: 'pkg-config', required: false, version: '>=4.2.0')

# preprocessor macros
add_project_arguments('-D CLAS_STRINGSPINNER_VERSION="' + meson.project_version() + '"', language: ['cpp'])
if hipo_dep.found()
  add_project_arguments('-D CLAS_STRINGSPINNER_HIPO', language: ['cpp'])
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
