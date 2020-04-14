macro(VerboseFindPackage package)
    if(${package}_ROOT)
        message(STATUS "--[${package}] Look for package in: ${${package}_ROOT}")
    else()
        message(STATUS "--[${package}] Look for package")
    endif()
    find_package(${ARGV})
    if(${package}_FOUND)
        message(STATUS "--[${package}] Package found.")
        if(${package}_VERSION)
            message(STATUS "--[${package}] Package version ${${package}_VERSION}")
        endif()
        if(${package}_INCLUDE_DIRS)
            message(STATUS "--[${package}] Package include dir: ${${package}_INCLUDE_DIRS}")
        endif()
        if(${package}_LIBRARY_DIRS)
            message(STATUS "--[${package}] Package lib dir: ${${package}_LIBRARY_DIRS}")
        endif()
        if(${package}_LIBRARIES)
            message(STATUS "--[${package}] Package libs: ${${package}_LIBRARIES}")
        endif()
    else()
        message(STATUS "--[${package}] Package NOT found")
    endif()
endmacro()
