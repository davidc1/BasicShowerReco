art_make(
          LIB_LIBRARIES
			uboone_BasicShowerReco_TwoDimTools
			art_Framework_Core
			art_Framework_Services_Registry
			art_Framework_Principal
			lardataobj_RecoBase
			${ART_FRAMEWORK_CORE}
			larcorealg_Geometry	
			lardata_Utilities # added now
		   	uboone_Geometry		
			larcore_Geometry_Geometry_service	
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
			${Boost_FILESYSTEM_LIBRARY}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_FRAMEWORK_SERVICES_OPTIONAL}
          TOOL_LIBRARIES  
			  uboone_BasicShowerReco_TwoDimTools
                          lardata_Utilities_LArFFT_service
                          larcorealg_Geometry
			  larevt_CalibrationDBI_IOVData
			  larevt_CalibrationDBI_Providers
			  lardataobj_RecoBase
			  larcore_Geometry_Geometry_service
			  lardata_Utilities
			  nutools_RandomUtils_NuRandomService_service
  			  ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			  art_Persistency_Common
			  art_Persistency_Provenance
			  ${ART_FRAMEWORK_SERVICES_BASIC}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			  art_Utilities
			  canvas
			  ${MF_MESSAGELOGGER}
			  ${MF_UTILITIES}
			  ${FHICLCPP}
			  ${CETLIB}
			  ${Boost_FILESYSTEM_LIBRARY}
			  ${Boost_SYSTEM_LIBRARY}
                          ${CLHEP}
			  ${ROOT_BASIC_LIB_LIST}
        )


install_headers()
install_fhicl()
install_source()

