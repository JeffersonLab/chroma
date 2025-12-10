function(CHECK_SANITIZER_OPTIONS SANITIZER_OPTION_LIST_IN OPTIONS_OUT)
 # Sanitizers only work for some compilers
 if( ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 3.1)
     OR ( CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 4.8)  )
    

  # If the value is enabled but no list is given, give both
 	if( "${SANITIZER_OPTION_LIST_IN}" STREQUAL "ON"
     OR "${SANITIZER_OPTION_LIST_IN}" STREQUAL "TRUE" )
		set(SANITIZER_OPT_STRING "-fsanitize=address,undefined")
  else()

		#if it is a list check each memeber
  	foreach(OPT IN ITEMS ${SANITIZER_OPTION_LIST_IN})
		  if ( ${OPT} STREQUAL "address"  ) 
		  	set(VALID TRUE)
	    elseif( ${OPT} STREQUAL "undefined" )
		  	set(VALID TRUE)
		  else()
		  	set(VALID FALSE)
		  endif()

      #If the member is valid add it on to the option
	  	if( VALID ) 
	  		if( NOT SANITIZER_OPT_STRING) 
	  			set(SANITIZER_OPT_STRING "-fsanitize=${OPT}")
				else()	
	  			set(SANITIZER_OPT_STRING ${SANITIZER_OPT_STRING},${OPT})
		    endif()
	    else() 
	      message(STATUS "Option ${OPT} is not a valid sanitizer option")
		  endif()
		endforeach()
   endif()
	else() 
    # if we don't know the compiler
	  message(STATUS "Don't know how to enable sanitizers for ${CMAKE_C_COMPILER}... Ignoring")
	  set(SANITIZER_OPT_STRING "")
	endif()
	set(${OPTIONS_OUT} ${SANITIZER_OPT_STRING} PARENT_SCOPE)
endfunction()
