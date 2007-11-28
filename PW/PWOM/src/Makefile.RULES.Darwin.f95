PhotoElectronPW_planet.o: PhotoElectronPW_planet.f
	${COMPILE.f77} ${Cflag1} -132 PhotoElectronPW_planet.f
IonHeatFluxPW_planet.o: IonHeatFluxPW_planet.f
	${COMPILE.f77} ${Cflag1} -132 IonHeatFluxPW_planet.f
collisionPW_planet.o: collisionPW_planet.f
	${COMPILE.f77} ${Cflag1} -132 collisionPW_planet.f
precipitationPW.o: precipitation_planet.f
	${COMPILE.f77} ${Cflag1} -132 precipitationPW_planet.f
glowex_planet.o: glowex_planet.f
	${COMPILE.f77} ${Cflag1} -132 glowex_planet.f
eHeatFluxPW.o: eHeatFluxPW.f
	${COMPILE.f77} ${Cflag1} -132 eHeatFluxPW.f
startupPW_planet.o: startupPW_planet.f
	${COMPILE.f77} ${Cflag1} -132 startupPW_planet.f
ssflux_planet.o: ssflux_planet.f
	${COMPILE.f77} ${Cflag1} -132 ssflux_planet.f
