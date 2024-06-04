project(
  'clas-stringspinner',
  'cpp',
  meson_version: '>=1.3',
  default_options: {
    'cpp_std':    'c++17',
    'buildtype':  'release',
    'libdir':     'lib',
  },
  version: '0.1.0',
)

# subprojects
stringspinner_proj = subproject('stringspinner')

# dependencies
fmt_dep           = dependency('fmt', method: 'pkg-config', version: '>= 9.1.0')
stringspinner_dep = dependency('stringspinner', fallback: ['stringspinner', 'stringspinner_dep'])

# config files
fs = import('fs')
etc_dir = get_option('sysconfdir') / meson.project_name()
add_project_arguments(
  '-DSTRINGSPINNER_ETCDIR="' + fs.relative_to(etc_dir, get_option('bindir')) + '"',
  '-DSTRINGSPINNER_PREFIX_ETCDIR="' + get_option('prefix') / etc_dir + '"',
  language: ['cpp'],
)
install_subdir(
  'config',
  install_dir: etc_dir,
  strip_directory: true,
  exclude_files: ['README.md'],
)

# main executable
main_exe = executable(
  'clas-stringspinner',
  'src/clas-stringspinner.cpp',
  dependencies: [ fmt_dep, stringspinner_dep ],
  install: true,
)

# test executable
test(
  'clas-stringspinner',
  main_exe,
  args: [
    '--num-events', '10000',
    '--pol-type', 'LU',
    '--beam-spin', 'p',
    '--cut-inclusive', '11,211,-211',
    '--seed', '828203788',
  ]
)

