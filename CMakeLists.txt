cmake_minimum_required(VERSION 3.0)

option( INSTALL_QPOSC_PLUGIN "Check to install qPOSC plugin" OFF )

# CloudCompare 'POSC' plugin
if (INSTALL_QPOSC_PLUGIN)
    project( QPOSC_PLUGIN )
    
    #load necessary libraries (see qPCV for an example)
    #add_subdirectory (LIB1)
    
    #if the plugin is an 'OpenGL filter', uncomment the line below
    #set( CC_OPENGL_FILTER ON BOOL)
    include( ../CMakePluginTpl.cmake )
    
    #set dependencies to necessary libraries (see qPCV for an example)
    #target_link_libraries( ${PROJECT_NAME} LIB1 )
    #include_directories( ${LIB1_INCLUDE_DIR} )
endif()
