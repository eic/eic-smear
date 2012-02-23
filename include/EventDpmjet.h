#ifndef _EIC_BUILDTREE_EVENTDPMJETH_
#define _EIC_BUILDTREE_EVENTDPMJETH_

//describles the event from DPMJET
#include <string>

#include <Rtypes.h>

#include "EventBase.h"

struct EventDpmjet : public EventBase {

	virtual bool Parse( const std::string& );

	Int_t ievent;
	Int_t I;
	Int_t process1;
	Int_t process2;
	Int_t IP;

	Int_t tgtparton;
	Int_t prjparton;
	Int_t nucleon;

	Double32_t xtgtparton;
	Double32_t xprjparton;
	Double32_t dtrueW2;
	Double32_t dtrueNu;
	Double32_t dtrueQ2;
	Double32_t dtrueY;
	Double32_t dtrueX;
	Double32_t theta_Evt;
	Double32_t photonFlux;
	ClassDef(EventDpmjet, 1)
};

#endif
