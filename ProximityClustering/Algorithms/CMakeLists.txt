art_make(
          LIB_LIBRARIES
		   	uboone_Geometry		
			lardata_Utilities	
			larcore_Geometry_Geometry_service	
			lardataobj_RecoBase
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_FRAMEWORK_SERVICES_OPTIONAL}

        )


install_headers()
install_fhicl()
install_source()

