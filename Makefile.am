# Options to pass to aclocal when aclocal.m4 is to be rebuilt by make
#ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libBuildTree.la

#include_HEADERS = $(top_srcdir)/include/*.h
include_HEADERS  =  \
$(top_srcdir)/include/Acceptance.h \
$(top_srcdir)/include/Armory.h \
$(top_srcdir)/include/BeamParticles.h \
$(top_srcdir)/include/Bremsstrahlung.h \
$(top_srcdir)/include/Detector.h \
$(top_srcdir)/include/Device.h \
$(top_srcdir)/include/DisKinematics.h \
$(top_srcdir)/include/EventBase.h \
$(top_srcdir)/include/EventDjangoh.h \
$(top_srcdir)/include/EventDpmjet.h \
$(top_srcdir)/include/EventFactory.h \
$(top_srcdir)/include/EventMCFilterABC.h \
$(top_srcdir)/include/EventMC.h \
$(top_srcdir)/include/EventMilou.h \
$(top_srcdir)/include/EventPepsi.h \
$(top_srcdir)/include/EventPythia.h \
$(top_srcdir)/include/EventRapgap.h \
$(top_srcdir)/include/EventS.h \
$(top_srcdir)/include/EventSmear.h \
$(top_srcdir)/include/File.h \
$(top_srcdir)/include/Forester.h \
$(top_srcdir)/include/functions.h \
$(top_srcdir)/include/Kinematics.h \
$(top_srcdir)/include/Particle.h \
$(top_srcdir)/include/ParticleIdentifier.h \
$(top_srcdir)/include/ParticleID.h \
$(top_srcdir)/include/ParticleMC.h \
$(top_srcdir)/include/ParticleMCSmeared.h \
$(top_srcdir)/include/Pid.h \
$(top_srcdir)/include/SmearEvent.h \
$(top_srcdir)/include/Smear.h \
$(top_srcdir)/include/VirtualEvent.h \
$(top_srcdir)/include/VirtualParticle.h

libBuildTree_la_SOURCES = \
$(top_srcdir)/src/Acceptance.cxx \
$(top_srcdir)/src/BeamParticles.cxx \
$(top_srcdir)/src/Bremsstrahlung.cxx \
$(top_srcdir)/src/BuildTree.cxx \
$(top_srcdir)/src/Detector.cxx \
$(top_srcdir)/src/Device.cxx \
$(top_srcdir)/src/EventDjangoh.cxx \
$(top_srcdir)/src/EventDpmjet.cxx \
$(top_srcdir)/src/EventFactory.cxx \
$(top_srcdir)/src/EventMC.cxx \
$(top_srcdir)/src/EventMilou.cxx \
$(top_srcdir)/src/EventPepsi.cxx \
$(top_srcdir)/src/EventPythia.cxx \
$(top_srcdir)/src/EventRapgap.cxx \
$(top_srcdir)/src/File.cxx \
$(top_srcdir)/src/Forester.cxx \
$(top_srcdir)/src/functions.cxx \
$(top_srcdir)/src/Kinematics.cxx \
$(top_srcdir)/src/ParticleIdentifier.cxx \
$(top_srcdir)/src/ParticleMC.cxx \
$(top_srcdir)/src/ParticleMCSmeared.cxx \
$(top_srcdir)/src/Pid.cxx \
$(top_srcdir)/src/Smear.cxx \
$(top_srcdir)/src/smear_event.cxx \
$(top_srcdir)/src/SmearTree.cxx \
$(top_srcdir)/cint/BuildTreeDict.cxx

libBuildTree_la_CPPFLAGS = \
-I$(top_srcdir)/include

libBuildTree_la_LIBADD = \
-lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -lm -ldl -lEG -lEGPythia6

libBuildTree_la_LDFLAGS = -version-info 2:0:0

if ENABLE_PYTHIA6_SUPPORT
include_HEADERS += \
$(top_srcdir)/include/Pythia6EventBuilder.h \
$(top_srcdir)/include/Pythia6.h \
$(top_srcdir)/include/Pythia6ParticleBuilder.h

libBuildTree_la_SOURCES += \
$(top_srcdir)/src/Pythia6.cxx \
$(top_srcdir)/src/Pythia6EventBuilder.cxx \
$(top_srcdir)/src/Pythia6ParticleBuilder.cxx

libBuildTree_la_LIBADD += @PYTHIA6_LIBS@

LINKDEF=$(top_srcdir)/cint/LinkDef.h
else
LINKDEF=$(top_srcdir)/cint/LinkDefNoPythia.h
endif

#$(top_srcdir)/cint/BuildTreeDict.cxx: \
#$(top_srcdir)/include/Acceptance.h \
#$(top_srcdir)/include/Armory.h \
#$(top_srcdir)/include/BeamParticles.h \
#$(top_srcdir)/include/Bremsstrahlung.h \
#$(top_srcdir)/include/Detector.h \
#$(top_srcdir)/include/Device.h \
#$(top_srcdir)/include/DisKinematics.h \
#$(top_srcdir)/include/EventBase.h \
#$(top_srcdir)/include/EventDjangoh.h \
#$(top_srcdir)/include/EventDpmjet.h \
#$(top_srcdir)/include/EventFactory.h \
#$(top_srcdir)/include/EventMCFilterABC.h \
#$(top_srcdir)/include/EventMC.h \
#$(top_srcdir)/include/EventMilou.h \
#$(top_srcdir)/include/EventPepsi.h \
#$(top_srcdir)/include/EventPythia.h \
#$(top_srcdir)/include/EventRapgap.h \
#$(top_srcdir)/include/EventS.h \
#$(top_srcdir)/include/EventSmear.h \
#$(top_srcdir)/include/File.h \
#$(top_srcdir)/include/Forester.h \
#$(top_srcdir)/include/functions.h \
#$(top_srcdir)/include/GeneratorType.h \
#$(top_srcdir)/include/Kinematics.h \
#$(top_srcdir)/include/Particle.h \
#$(top_srcdir)/include/ParticleIdentifier.h \
#$(top_srcdir)/include/ParticleID.h \
#$(top_srcdir)/include/ParticleMC.h \
#$(top_srcdir)/include/ParticleMCSmeared.h \
#$(top_srcdir)/include/Pid.h \
#$(top_srcdir)/include/SmearEvent.h \
#$(top_srcdir)/include/Smear.h \
#$(top_srcdir)/include/VirtualEvent.h \
#$(top_srcdir)/include/VirtualParticle.h \
#$(top_srcdir)/cint/LinkDef.h
#	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

if ENABLE_PYTHIA6_SUPPORT
test:
	@echo "User enabled PYTHIA 6"
#	@echo @enable_pythia6@
	@echo @PYTHIA6_LIBS@
	@echo $(libBuildTree_la_LIBADD)
else
test:
	@echo Not using PYTHIA 6
	@echo $(libBuildTree_la_LIBADD)
endif

$(top_srcdir)/cint/BuildTreeDict.cxx: $(include_HEADERS) $(LINKDEF)
	@echo selected linkdef $^
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f $(BUILT_SOURCES) cint/*Dict*


