
add_library(virtualAmps   "bsyspinors.f"
                          "bsyspinit.f"
                          "ampggint.f"
                          "ampgg.f"
                          "ampqqint.f"
                          "ampqq.f"
                          "bsypptt.f")


add_dependencies(virtualAmps ff)
add_dependencies(virtualAmps qcdloop)
target_link_libraries(virtualAmps "-lqcdloop" "-lff")