function (print VAR)
  set(NAME ${VAR})
  list (APPEND VARS ${${VAR}})

  set(MSG)
  set(IS_FIRST TRUE)

  foreach(ELM ${VARS})
    if (IS_FIRST)
      set(MSG ${ELM})
      set(IS_FIRST FALSE)
    else ()
      set(MSG ${MSG}, ${ELM})
    endif ()
  endforeach()

  message("${NAME} = "${MSG})
endfunction(print VAR)

