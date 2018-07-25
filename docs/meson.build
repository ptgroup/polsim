config = configuration_data()
config.set('PROJECT_NAME', meson.project_name())
config.set('VERSION', meson.project_version())
config.set('SRCDIR', meson.source_root())
config.set('DESTDIR', join_paths(meson.build_root(), 'docs'))
config.set('README_PATH', join_paths(meson.source_root(), 'README.md'))

doxyfile = configure_file(
    input : 'Doxyfile.in',
    output : 'Doxyfile',
    configuration : config
)

custom_target(
    'docs',
    input : doxyfile,
    output : 'docs',
    command : ['doxygen', '@INPUT@'],
    build_by_default : true
)
