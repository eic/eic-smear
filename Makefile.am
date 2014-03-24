lib_LTLIBRARIES = libeicsmear.la

# Note we list headers separately for each subdirectory so
# we can match this structure when doing make install.

# Headers in eicsmear/erhic.
erhic_headers = \
include/eicsmear/functions.h \
include/eicsmear/erhic/BeamParticles.h \
include/eicsmear/erhic/DisKinematics.h \
include/eicsmear/erhic/EventBase.h \
include/eicsmear/erhic/EventDis.h \
include/eicsmear/erhic/EventDjangoh.h \
include/eicsmear/erhic/EventDpmjet.h \
include/eicsmear/erhic/EventFactory.h \
include/eicsmear/erhic/EventGmcTrans.h \
include/eicsmear/erhic/EventMCFilterABC.h \
include/eicsmear/erhic/EventMC.h \
include/eicsmear/erhic/EventMilou.h \
include/eicsmear/erhic/EventPepsi.h \
include/eicsmear/erhic/EventPythia.h \
include/eicsmear/erhic/EventRapgap.h \
include/eicsmear/erhic/File.h \
include/eicsmear/erhic/Forester.h \
include/eicsmear/erhic/Kinematics.h \
include/eicsmear/erhic/Particle.h \
include/eicsmear/erhic/ParticleIdentifier.h \
include/eicsmear/erhic/ParticleMC.h \
include/eicsmear/erhic/Pid.h \
include/eicsmear/erhic/VirtualEvent.h \
include/eicsmear/erhic/VirtualParticle.h

# Headers in eicsmear/smear.
smear_headers = \
include/eicsmear/smear/Acceptance.h \
include/eicsmear/smear/Bremsstrahlung.h \
include/eicsmear/smear/Detector.h \
include/eicsmear/smear/Device.h \
include/eicsmear/smear/Distributor.h \
include/eicsmear/smear/EventDisFactory.h \
include/eicsmear/smear/EventFactory.h \
include/eicsmear/smear/EventS.h \
include/eicsmear/smear/EventSmear.h \
include/eicsmear/smear/FormulaString.h \
include/eicsmear/smear/functions.h \
include/eicsmear/smear/ParticleID.h \
include/eicsmear/smear/ParticleMCS.h \
include/eicsmear/smear/PerfectID.h \
include/eicsmear/smear/PlanarTracker.h \
include/eicsmear/smear/RadialTracker.h \
include/eicsmear/smear/Smear.h \
include/eicsmear/smear/Smearer.h \
include/eicsmear/smear/Tracker.h

# Headers in eicsmear/hadronic (PYTHIA dependent).
hadronic_headers = \
include/eicsmear/hadronic/EventMC.h \
include/eicsmear/hadronic/EventPythia.h \
include/eicsmear/hadronic/ParticleMC.h \
include/eicsmear/hadronic/Pythia6EventFactory.h \
include/eicsmear/hadronic/EventSmear.h

# PYTHIA-dependent headers also in eicsmear/erhic.
pythia_headers = \
include/eicsmear/erhic/Pythia6EventBuilder.h \
include/eicsmear/erhic/Pythia6.h \
include/eicsmear/erhic/Pythia6ParticleBuilder.h

# All non-PYTHIA-dependent source files.
libeicsmear_la_SOURCES = \
src/erhic/BeamParticles.cxx \
src/erhic/BuildTree.cxx \
src/erhic/EventDis.cxx \
src/erhic/EventDjangoh.cxx \
src/erhic/EventDpmjet.cxx \
src/erhic/EventFactory.cxx \
src/erhic/EventGmcTrans.cxx \
src/erhic/EventMC.cxx \
src/erhic/EventMilou.cxx \
src/erhic/EventPepsi.cxx \
src/erhic/EventPythia.cxx \
src/erhic/EventRapgap.cxx \
src/erhic/File.cxx \
src/erhic/Forester.cxx \
src/erhic/Kinematics.cxx \
src/erhic/ParticleIdentifier.cxx \
src/erhic/ParticleMC.cxx \
src/erhic/Pid.cxx \
src/functions.cxx \
src/smear/Acceptance.cxx \
src/smear/Bremsstrahlung.cxx \
src/smear/Detector.cxx \
src/smear/Device.cxx \
src/smear/Distributor.cxx \
src/smear/EventDisFactory.cxx \
src/smear/EventSmear.cxx \
src/smear/FormulaString.cxx \
src/smear/ParticleID.cxx \
src/smear/ParticleMCS.cxx \
src/smear/PerfectID.cxx \
src/smear/PlanarTracker.cxx \
src/smear/RadialTracker.cxx \
src/smear/Smear.cxx \
src/smear/SmearTree.cxx \
src/smear/Tracker.cxx \
cint/LinkDef.h \
cint/BuildTreeDict.cxx \
cint/smearLinkDef.h \
cint/smearDict.cxx

# Use this construction so when we call make install
# all the files end up in the correct subdirectories.
eicsmearincludedir = $(includedir)/eicsmear
eicsmearinclude_HEADERS = include/eicsmear/functions.h
erhicincludedir = $(includedir)/eicsmear/erhic
erhicinclude_HEADERS = $(erhic_headers)
smearincludedir = $(includedir)/eicsmear/smear
smearinclude_HEADERS = $(smear_headers)

# Set default preprocessor flags.
libeicsmear_la_CPPFLAGS = -I$(top_srcdir)/include

# Link various ROOT libraries.
libeicsmear_la_LIBADD = \
-lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
-lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -lm -ldl -lEG

# ------------------------------------------------
# If the user passed a PYTHIA6 library directory,
# build the PYTHIA6-dependent code.
if ENABLE_PYTHIA6_SUPPORT
# Add directories for PYTHIA-dependent headers (both
# in eicsmear/hadronic and eicsmear/erhic).
hadronicincludedir = $(includedir)/eicsmear/hadronic
hadronicinclude_HEADERS = $(hadronic_headers)
pythiaincludedir = $(includedir)/eicsmear/erhic
pythiainclude_HEADERS = $(pythia_headers)

# Add the PYTHIA-dependent source files.
libeicsmear_la_SOURCES += \
src/erhic/Pythia6.cxx \
src/erhic/Pythia6EventBuilder.cxx \
src/erhic/Pythia6ParticleBuilder.cxx \
src/hadronic/EventMC.cxx \
src/hadronic/EventPythia.cxx \
src/hadronic/EventSmear.cxx \
src/hadronic/ParticleMC.cxx \
src/hadronic/Pythia6EventFactory.cxx \
cint/hadronicLinkDef.h \
cint/hadronicDict.cxx \
cint/LinkDefPythia.h \
cint/pythiaDict.cxx


# In the smearing routine we need to be able to turn
# on and off support for smearing hadronic events.
# Turn on WITH_PYTHIA6 to allow conditional compilation.
libeicsmear_la_CPPFLAGS += -DWITH_PYTHIA6

# Link against the PYTHIA library.
libeicsmear_la_LIBADD += @PYTHIA6_LIBS@
endif # ENABLE_PYTHIA6_SUPPORT
# ------------------------------------------------

# Build ROOT CINT dictionaries separately for each component.
cint/hadronicDict.cxx: $(hadronic_headers) cint/hadronicLinkDef.h
	rootcint -f $@ -c -I$(top_srcdir)/include $^

cint/pythiaDict.cxx: $(pythia_headers) cint/LinkDefPythia.h
	rootcint -f $@ -c -I$(top_srcdir)/include $^

cint/smearDict.cxx: $(smear_headers) cint/smearLinkDef.h
	rootcint -f $@ -c -I$(top_srcdir)/include $^

cint/BuildTreeDict.cxx: $(erhic_headers) cint/LinkDef.h
	rootcint -f $@ -c -I$(top_srcdir)/include $^

EXTRA_DIST = CMakeLists.txt cmake
dist-hook:
	rm -rf `find $(distdir) -name "*Dict*"`

clean-local:
	rm -f $(BUILT_SOURCES) cint/*Dict* .libs _libs
