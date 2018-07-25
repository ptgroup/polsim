project('polsim', 'cpp', default_options : ['cpp_std=c++14', 'warnlevel=3'])

boost = dependency('boost')
incdir = include_directories('include')
subdir('src')

doxygen = find_program('doxygen', required : false)
if doxygen.found()
    message('Found doxygen; adding documentation target.')
    subdir('docs')
else
    message('Doxygen not installed; no documentation will be generated.')
endif
