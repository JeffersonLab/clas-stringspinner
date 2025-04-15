project(
  'clas-stringspinner',
  'cpp',
  meson_version: '>=1.3',
  default_options: {
    'cpp_std':    'c++17',
    'buildtype':  'release',
    'libdir':     'lib',
    'wrap_mode':  'forcefallback',
  },
  version: '0.1.4',
)

# subprojects
stringspinner_proj = subproject('stringspinner')

# dependencies
fmt_dep           = dependency('fmt', method: 'pkg-config', version: '>= 9.1.0')
stringspinner_dep = dependency('stringspinner', fallback: ['stringspinner', 'stringspinner_dep'])
pythia_dep        = dependency('pythia', fallback: ['stringspinner', 'stringspinner_pythia_dep'], version: '8.312')

# main executable
add_project_arguments('-DCLAS_STRINGSPINNER_VERSION="' + meson.project_version() + '"', language: ['cpp'])
main_exe = executable(
  meson.project_name(),
  'src' / meson.project_name() + '.cpp',
  'src/config/src/common.cc',
  include_directories: include_directories('src/config'),
  dependencies: [ fmt_dep, stringspinner_dep, pythia_dep ],
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

