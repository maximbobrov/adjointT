TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11

QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11
