#1)Primero se incluyen los codigos

	#Incluye el directorio con tus propios codigos
	INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

	#Incluye completamente la libreria ParadisEO
	INCLUDE_DIRECTORIES(${PARADISEO_INCLUDE_DIR})			
				    
#2)Define los ejecutables y dependencias
	
	ADD_EXECUTABLE(AGS-permuta AGS-permuta.cpp)
	ADD_DEPENDENCIES(AGS-permuta lDiscovery)
	ADD_EXECUTABLE(AGS-island AGS-island.cpp)
	ADD_DEPENDENCIES(AGS-island lDiscovery)

#3)Vincular las librerias para el ejecutable

	TARGET_LINK_LIBRARIES(AGS-permuta ${PARADISEO_LIBRARIES} lDiscovery)
	TARGET_LINK_LIBRARIES(AGS-island ${PARADISEO_LIBRARIES} lDiscovery)
