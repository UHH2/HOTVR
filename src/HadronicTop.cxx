#include "UHH2/HOTVR/include/HadronicTop.h"

using namespace uhh2;
using namespace std;

HadronicTop::HadronicTop(uhh2::Context &ctx, const std::string & handle_name) :
  h_tophad(ctx.get_handle<vector<GenTop> >(handle_name)){}

bool HadronicTop::process(Event &event) {
  vector<GenTop> gentops;
  const vector<GenParticle> genparticles = *event.genparticles;
  for (const GenParticle & genp : genparticles)
    {
      if (abs(genp.pdgId()) == 6)
	{
	  GenTop top;
	  const GenParticle *w = genp.daughter(&genparticles, 1);
	  const GenParticle *b = genp.daughter(&genparticles, 2);
	  if (!w || !b){
	    //if(throw_on_failure) throw runtime_error("TTbarGen: top has not ==2 daughters");
	    continue;
	  }
	  if (abs(w->pdgId()) != 24){
	    std::swap(w, b);
	  }
	  if (abs(w->pdgId()) != 24){
	    //if(throw_on_failure) throw runtime_error("TTbarGen: top has no W daughter");
	    continue;
	  }
	  const GenParticle *wd1 = w->daughter(&genparticles, 1);
	  const GenParticle *wd2 = w->daughter(&genparticles, 2);
	  if ( !( (abs(wd1->pdgId()) < 9) && (abs(wd2->pdgId()) < 9) ) )continue;
	  top.set_top(genp.v4());
	  top.set_b(b->v4());
	  top.set_W(w->v4());
	  top.set_q1(wd1->v4());
	  top.set_q2(wd2->v4());
	  gentops.push_back(top);
	}
    }
  event.set(h_tophad, gentops);
  return true;
}

HadronicTop::~HadronicTop() {}
