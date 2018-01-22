#1)Primero se incluyen los directorios que uno necesita

	#Incluye completamente la libreria ParadisEO
	INCLUDE_DIRECTORIES(${PARADISEO_INCLUDE_DIR})

	#Incluye el directorio con tus propios codigos
	INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

#2)Se debe establecer que la libreria debe ser compilada
	
	SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

#3)Se deben agregar los propios codigos ".cpp" (uno elige el nombre de la lista)

	SET(CODIGOS
			localizacionEvalPenal.cpp
			escenario.cpp
	)

#4)Declara tu propia libreria: todos los archivos ".cpp" que estan definidos en la lista (uno puede elegir el nombre de la libreria)

	ADD_LIBRARY(lDiscovery STATIC ${CODIGOS}) 

