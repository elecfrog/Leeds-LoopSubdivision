cmake_minimum_required(VERSION 3.25)
project(Leeds-LoopSubdivision)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

set(CMAKE_PREFIX_PATH "C://Qt//6.5.1//msvc2019_64")

find_package(Qt6 COMPONENTS
        Core
        Gui
        Widgets
        OpenGLWidgets
        REQUIRED)

find_package(OpenGL REQUIRED)


add_executable(Leeds-LoopSubdivision
        src/ArcBall.h
        src/ArcBallWidget.h
        src/Cartesian3.h
        src/DirectedEdgeSurface.h
        src/Homogeneous4.h
        src/Matrix4.h
        src/Quaternion.h
        src/RGBAImage.h
        src/RGBAValue.h
        src/RenderController.h
        src/RenderParameters.h
        src/RenderWidget.h
        src/RenderWindow.h
        src/ArcBall.cpp
        src/ArcBallWidget.cpp
        src/Cartesian3.cpp
        src/DirectedEdgeSurface.cpp
        src/Homogeneous4.cpp
        src/Matrix4.cpp
        src/Quaternion.cpp
        src/RGBAImage.cpp
        src/RGBAValue.cpp
        src/RenderController.cpp
        src/RenderWidget.cpp
        src/RenderWindow.cpp
        src/SphereVertices.cpp
        src/main.cpp
)

target_link_libraries(Leeds-LoopSubdivision
        PRIVATE
        Qt::Core
        Qt::Gui
        Qt::Widgets
        Qt::OpenGLWidgets

        OpenGL::GL
)

if (WIN32 AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(DEBUG_SUFFIX)
    if (MSVC AND CMAKE_BUILD_TYPE MATCHES "Debug")
        set(DEBUG_SUFFIX "d")
    endif ()
    set(QT_INSTALL_PATH "${CMAKE_PREFIX_PATH}")
    if (NOT EXISTS "${QT_INSTALL_PATH}/bin")
        set(QT_INSTALL_PATH "${QT_INSTALL_PATH}/..")
        if (NOT EXISTS "${QT_INSTALL_PATH}/bin")
            set(QT_INSTALL_PATH "${QT_INSTALL_PATH}/..")
        endif ()
    endif ()
    if (EXISTS "${QT_INSTALL_PATH}/plugins/platforms/qwindows${DEBUG_SUFFIX}.dll")
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory
                "$<TARGET_FILE_DIR:${PROJECT_NAME}>/plugins/platforms/")
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                "${QT_INSTALL_PATH}/plugins/platforms/qwindows${DEBUG_SUFFIX}.dll"
                "$<TARGET_FILE_DIR:${PROJECT_NAME}>/plugins/platforms/")
    endif ()
    foreach (QT_LIB Core Gui Widgets Sql)
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                "${QT_INSTALL_PATH}/bin/Qt6${QT_LIB}${DEBUG_SUFFIX}.dll"
                "$<TARGET_FILE_DIR:${PROJECT_NAME}>")
    endforeach (QT_LIB)
endif ()
