#pragma once

#include "Lund.h"

#include <hipo4/writer.h>

namespace css {

  class Hipo {

    private:
#ifdef CLAS_STRINGSPINNER_HIPO
      hipo::writer writer;
      hipo::event event;
      hipo::schema s_mc_header{"MC::Header", 40, 0};
      hipo::schema s_mc_event{"MC::Event", 40, 1};
      hipo::schema s_mc_particle{"MC::Particle", 40, 2};
      hipo::schema s_mc_lund{"MC::Lund", 40, 3};
      hipo::bank b_mc_header;
      hipo::bank b_mc_event;
      hipo::bank b_mc_particle;
      hipo::bank b_mc_lund;
#endif

    public: // ==================================================================================

      Hipo(std::string hipo_file_name)
      {
#ifdef CLAS_STRINGSPINNER_HIPO
        s_mc_header.parse("run/I,event/I,type/B,helicity/F");
        s_mc_event.parse("npart/S,atarget/S,ztarget/S,ptarget/F,pbeam/F,btype/S,ebeam/F,targetid/S,processid/S,weight/F");
        s_mc_particle.parse("pid/I,px/F,py/F,pz/F,vx/F,vy/F,vz/F,vt/F");
        s_mc_lund.parse("index/B,lifetime/F,type/B,pid/I,parent/B,daughter/B,px/F,py/F,pz/F,energy/F,mass/F,vx/F,vy/F,vz/F");
        writer.getDictionary().addSchema(s_mc_header);
        writer.getDictionary().addSchema(s_mc_event);
        writer.getDictionary().addSchema(s_mc_particle);
        writer.getDictionary().addSchema(s_mc_lund);
        b_mc_header   = hipo::bank(s_mc_header);
        b_mc_event    = hipo::bank(s_mc_event);
        b_mc_particle = hipo::bank(s_mc_particle);
        b_mc_lund     = hipo::bank(s_mc_lund);
        writer.open(hipo_file_name.c_str());
#else
        throw std::runtime_error("HIPO file output was not enabled for this build of clas-stringspinner");
#endif
      }

      // ==================================================================================

      void Stream(
          LundHeader const& head,
          std::vector<LundParticle> const& pars,
          css::evnum_t evnum
          )
      {
#ifdef CLAS_STRINGSPINNER_HIPO
        // fill `MC::Header`
        b_mc_header.setRows(1);
        b_mc_header.putInt("run",        0, 11);
        b_mc_header.putInt("event",      0, static_cast<int>(evnum));
        b_mc_header.putByte("type",      0, 0);
        b_mc_header.putFloat("helicity", 0, head.beam_spin);
        // fill `MC::Event`
        b_mc_event.setRows(1);
        b_mc_event.putShort("npart",     0, head.num_particles);
        b_mc_event.putShort("atarget",   0, head.target_mass);
        b_mc_event.putShort("ztarget",   0, head.target_atomic_num);
        b_mc_event.putFloat("ptarget",   0, head.target_spin);
        b_mc_event.putFloat("pbeam",     0, head.beam_spin);
        b_mc_event.putShort("btype",     0, head.beam_type);
        b_mc_event.putFloat("ebeam",     0, head.beam_energy);
        b_mc_event.putShort("targetid",  0, head.nucleon_pdg);
        b_mc_event.putShort("processid", 0, head.process_id);
        b_mc_event.putFloat("weight",    0, head.event_weight);
        // fill `MC::Particle`
        std::vector<LundParticle> pars_final; // `pars` with `Pythia::Particle::isFinal()` (i.e. `type==1`)
        for(auto par : pars)
          if(par.type == 1)
            pars_final.push_back(par);
        b_mc_particle.setRows(pars_final.size());
        int i_par = 0;
        for(auto const& par : pars_final) {
          b_mc_particle.putInt("pid",  i_par, par.pdg);
          b_mc_particle.putFloat("px", i_par, par.px);
          b_mc_particle.putFloat("py", i_par, par.py);
          b_mc_particle.putFloat("pz", i_par, par.pz);
          b_mc_particle.putFloat("vx", i_par, par.vx);
          b_mc_particle.putFloat("vy", i_par, par.vy);
          b_mc_particle.putFloat("vz", i_par, par.vz);
          b_mc_particle.putFloat("vt", i_par, 0);
          i_par++;
        }
        // fill `MC::Lund`
        b_mc_lund.setRows(pars.size());
        i_par = 0;
        for(auto const& par : pars) {
          b_mc_lund.putByte("index",     i_par, par.index);
          b_mc_lund.putFloat("lifetime", i_par, par.lifetime);
          b_mc_lund.putByte("type",      i_par, par.type);
          b_mc_lund.putInt("pid",        i_par, par.pdg);
          b_mc_lund.putByte("parent",    i_par, par.mother1);
          b_mc_lund.putByte("daughter",  i_par, par.daughter1);
          b_mc_lund.putFloat("px",       i_par, par.px);
          b_mc_lund.putFloat("py",       i_par, par.py);
          b_mc_lund.putFloat("pz",       i_par, par.pz);
          b_mc_lund.putFloat("energy",   i_par, par.energy);
          b_mc_lund.putFloat("mass",     i_par, par.mass);
          b_mc_lund.putFloat("vx",       i_par, par.vx);
          b_mc_lund.putFloat("vy",       i_par, par.vy);
          b_mc_lund.putFloat("vz",       i_par, par.vz);
          i_par++;
        }
        // write
        event.reset();
        event.addStructure(b_mc_header);
        event.addStructure(b_mc_event);
        event.addStructure(b_mc_particle);
        event.addStructure(b_mc_lund);
        writer.addEvent(event);
        // reset
        b_mc_header.reset();
        b_mc_event.reset();
        b_mc_particle.reset();
        b_mc_lund.reset();
#endif
      }

      // ==================================================================================

      void Close()
      {
#ifdef CLAS_STRINGSPINNER_HIPO
        writer.close();
#endif
      }

  };
}
