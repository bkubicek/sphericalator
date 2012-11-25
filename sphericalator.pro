#-------------------------------------------------
#
# Project created by QtCreator 2012-11-25T08:42:53
#
#-------------------------------------------------

QT       += core

QT       += gui opengl

TARGET = sphericalator
#CONFIG   += console
CONFIG   -= app_bundle
QMAKE_CXXFLAGS+= -fopenmp -O3 -march=pentium4   -ffast-math
QMAKE_LFLAGS += -openmp
LIBS += -fopenmp -fomit-frame-pointer

TEMPLATE = app


SOURCES += main.cpp \
    window.cpp \
    glwidget.cpp \
    sphericalcalculator.cpp \
    patch.cpp

HEADERS += \
    window.h \
    glwidget.h \
    sphericalcalculator.h \
    patch.h










