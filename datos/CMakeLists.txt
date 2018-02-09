######################################################################################
### Copia los archivos "param" a la carpeta de compilación para un facil uso
### 	UNIX
### --> ejecutar "make instalar" para copiar los archivos en el directorio donde se
###	construye la aplicación
###
###	Windows
### --> ejecutar "mingw32-make instalar" para copiar los archivos en el directorio donde se
###	construye la aplicación
######################################################################################

add_custom_target(instalar DEPENDS 
			${CMAKE_SOURCE_DIR}/datos/AGS-permuta.param
			${CMAKE_SOURCE_DIR}/datos/Anclas1.txt
			${CMAKE_SOURCE_DIR}/datos/Anclas2.txt
			${CMAKE_SOURCE_DIR}/datos/Anclas3.txt
			${CMAKE_SOURCE_DIR}/datos/Allposition1.txt
			${CMAKE_SOURCE_DIR}/datos/Allposition2.txt
			${CMAKE_SOURCE_DIR}/datos/Allposition3.txt
			${CMAKE_SOURCE_DIR}/datos/Matriz1.txt
			${CMAKE_SOURCE_DIR}/datos/Matriz2.txt
			${CMAKE_SOURCE_DIR}/datos/Matriz3.txt
				  ${CMAKE_SOURCE_DIR}/datos)

add_custom_command(
    TARGET instalar
    POST_BUILD
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/AGS-permuta.param
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Anclas1.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Anclas2.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Anclas3.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Allposition1.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Allposition2.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Allposition3.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Matriz1.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Matriz2.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

    ARGS -E copy_if_different
      ${CMAKE_SOURCE_DIR}/datos/Matriz3.txt
      ${CMAKE_BINARY_DIR}/application
    COMMAND ${CMAKE_COMMAND}

)

