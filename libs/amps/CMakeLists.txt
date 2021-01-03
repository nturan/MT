add_library(amps "gg_ttb.f"
                 "gg_ttbg.f"
                 "gu_ttbu.f"
                 "gub_ttbub.f"
                 "ubg_ttbub.f"
                 "ubu_ttb.f"
                 "ubu_ttbg.f"
                 "ug_ttbu.f"
                 "uub_ttb.f"
                 "uub_ttbg.f")

target_link_libraries(amps "-ldhelas")