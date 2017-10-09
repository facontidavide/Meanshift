cmake_minimum_required(VERSION 2.8)

project(MeanShiftDemo)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
SET(CMAKE_AUTOUIC ON)

find_package(Qt5 COMPONENTS Core Gui Widgets PrintSupport)

include_directories(
    include
    demo
    ${Qt5Core_INCLUDE_DIRS}
    ${Qt5Gui_INCLUDE_DIRS}
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5PrintSupport_INCLUDE_DIRS}
    )

add_definitions( -std=c++11 -fPIC -march=native)

add_executable(${PROJECT_NAME}
    demo/main.cpp
    include/MeanShift/MeanShift.h

    demo/mainwindow.h
    demo/mainwindow.cpp
    demo/mainwindow.ui

    demo/qcustomplot.h
    demo/qcustomplot.cpp
    )

target_link_libraries(${PROJECT_NAME}
    Qt5::Core
    Qt5::Widgets
    Qt5::PrintSupport
)
