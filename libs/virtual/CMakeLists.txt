
add_library(virtualAmps   "bsyspinors.f"
                          "bsyspinit.f"
                          "ampggint.f"
                          "ampgg.f"
                          "ampqqint.f"
                          "ampqq.f"
                          "bsypptt.f")

target_link_libraries(virtualAmps "-lqcdloop" "-lff")