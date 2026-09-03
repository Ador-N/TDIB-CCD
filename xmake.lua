set_project("subdiv-ccd")
set_languages("cxxlatest")
set_optimize("fastest")

add_rules("mode.debug", "mode.release")
add_requires("eigen","fmt")

target("scene")
    set_kind("binary")
    add_includedirs("core", "external/openGJK/scalar/include", {public = true})
    add_packages("eigen","fmt", {public = true})
    add_files("scene/demo.cpp","core/*.cpp", "external/openGJK/scalar/openGJK.c")
