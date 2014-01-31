######################################################################
# Automatically generated by qmake (2.01a) pt. sty 31 13:09:40 2014
######################################################################

TEMPLATE = app
TARGET = 
DEPENDPATH += . engine sim
INCLUDEPATH += . engine sim
QT += opengl

# Input
HEADERS += engine/dpintegrator.h \
           engine/entity.h \
           engine/geometry.h \
           engine/kerr.h \
           engine/kerr_coords.h \
           engine/numeric.h \
           engine/particle.h \
           engine/rk4integrator.h \
           engine/schw.h \
           sim/interface.h \
           sim/simulation.h
SOURCES += engine/dpintegrator.cpp \
           engine/entity.cpp \
           engine/geometry.cpp \
           engine/kerr.cpp \
           engine/kerr_coords.cpp \
           engine/numeric.cpp \
           engine/particle.cpp \
           engine/rk4integrator.cpp \
           engine/schw.cpp \
           sim/interface.cpp \
           sim/main.cpp \
           sim/simulation.cpp
