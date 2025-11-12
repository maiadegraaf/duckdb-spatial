add_subdirectory(util)
add_subdirectory(modules)
add_subdirectory(geometry)
add_subdirectory(index)
add_subdirectory(operators)

set(EXTENSION_SOURCES
    ${EXTENSION_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/spatial_extension.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_types.cpp
    PARENT_SCOPE)
