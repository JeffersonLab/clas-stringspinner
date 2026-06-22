#pragma once

#include "Lund.h"

#ifdef CLAS_STRINGSPINNER_HIPO
#include <hipo4/writer.h>
#endif

namespace string_spinner {

  class Hipo {

    int const MC_RUN_NUM = 11;

    private:
#ifdef CLAS_STRINGSPINNER_HIPO
      hipo::writer writer;
      hipo::event event;
      hipo::schema s_mc_header{"MC::Header", 40, 0};
      hipo::schema s_mc_event{"MC::Event", 40, 1};
      hipo::schema s_mc_particle{"MC::Particle", 40, 2};
      hipo::schema s_mc_lund{"MC::Lund", 40, 3};
      hipo::schema s_mc_user{"MC::User", 40, 5};
      hipo::bank b_mc_header;
      hipo::bank b_mc_event;
      hipo::bank b_mc_particle;
      hipo::bank b_mc_lund;
      hipo::bank b_mc_user;
#endif

    public: // ==================================================================================

      Hipo(std::string hipo_file_name)
      {
#ifdef CLAS_STRINGSPINNER_HIPO
        s_mc_header.parse("run/I,event/I,type/B,helicity/F");
        s_mc_event.parse("npart/S,atarget/S,ztarget/S,ptarget/F,pbeam/F,btype/S,ebeam/F,targetid/S,processid/S,weight/F");
        s_mc_particle.parse("pid/I,px/F,py/F,pz/F,vx/F,vy/F,vz/F,vt/F");
        s_mc_lund.parse("index/B,lifetime/F,type/B,pid/I,parent/B,daughter/B,px/F,py/F,pz/F,energy/F,mass/F,vx/F,vy/F,vz/F");
        s_mc_user.parse("userVar/F");
        writer.getDictionary().addSchema(s_mc_header);
        writer.getDictionary().addSchema(s_mc_event);
        writer.getDictionary().addSchema(s_mc_particle);
        writer.getDictionary().addSchema(s_mc_lund);
        writer.getDictionary().addSchema(s_mc_user);
        b_mc_header   = hipo::bank(s_mc_header);
        b_mc_event    = hipo::bank(s_mc_event);
        b_mc_particle = hipo::bank(s_mc_particle);
        b_mc_lund     = hipo::bank(s_mc_lund);
        b_mc_user     = hipo::bank(s_mc_user);
        writer.open(hipo_file_name.c_str());
#else
        throw std::runtime_error("HIPO file output was not enabled for this build of clas-stringspinner");
#endif
      }

      // ==================================================================================

      void Stream(
          LundHeader const& head,
          std::vector<LundParticle> const& pars,
          string_spinner::evnum_t evnum
          )
      {
#ifdef CLAS_STRINGSPINNER_HIPO
        // fill `MC::Header`
        b_mc_header.setRows(1);
        b_mc_header.putInt("run",        0, MC_RUN_NUM);
        b_mc_header.putInt("event",      0, static_cast<int32_t>(evnum));
        b_mc_header.putByte("type",      0, 0);
        b_mc_header.putFloat("helicity", 0, static_cast<float>(head.beam_spin));
        // fill `MC::Event`
        b_mc_event.setRows(1);
        b_mc_event.putShort("npart",     0, static_cast<int16_t>(head.num_particles));
        b_mc_event.putShort("atarget",   0, static_cast<int16_t>(head.target_mass_num));
        b_mc_event.putShort("ztarget",   0, static_cast<int16_t>(head.target_atomic_num));
        b_mc_event.putFloat("ptarget",   0, static_cast<float>(head.target_spin));
        b_mc_event.putFloat("pbeam",     0, static_cast<float>(head.beam_spin));
        b_mc_event.putShort("btype",     0, static_cast<int16_t>(head.beam_type));
        b_mc_event.putFloat("ebeam",     0, static_cast<float>(head.beam_energy));
        b_mc_event.putShort("targetid",  0, static_cast<int16_t>(head.nucleon_pdg));
        b_mc_event.putShort("processid", 0, static_cast<int16_t>(head.process_id));
        b_mc_event.putFloat("weight",    0, static_cast<float>(head.event_weight));
        // fill `MC::Particle`
        std::vector<LundParticle> pars_final; // `pars` with `Pythia::Particle::isFinal()` (i.e. `type==1`)
        for(auto par : pars)
          if(par.type == 1)
            pars_final.push_back(par);
        b_mc_particle.setRows(pars_final.size());
        int i_par = 0;
        for(auto const& par : pars_final) {
          b_mc_particle.putInt("pid",  i_par, static_cast<int32_t>(par.pdg));
          b_mc_particle.putFloat("px", i_par, static_cast<float>(par.px));
          b_mc_particle.putFloat("py", i_par, static_cast<float>(par.py));
          b_mc_particle.putFloat("pz", i_par, static_cast<float>(par.pz));
          b_mc_particle.putFloat("vx", i_par, static_cast<float>(par.vx));
          b_mc_particle.putFloat("vy", i_par, static_cast<float>(par.vy));
          b_mc_particle.putFloat("vz", i_par, static_cast<float>(par.vz));
          b_mc_particle.putFloat("vt", i_par, 0);
          i_par++;
        }
        // fill `MC::Lund`
        b_mc_lund.setRows(pars.size());
        i_par = 0;
        for(auto const& par : pars) {
          b_mc_lund.putByte("index",     i_par, static_cast<int8_t>(par.index));
          b_mc_lund.putFloat("lifetime", i_par, static_cast<float>(par.lifetime));
          b_mc_lund.putByte("type",      i_par, static_cast<int8_t>(par.type));
          b_mc_lund.putInt("pid",        i_par, static_cast<int32_t>(par.pdg));
          b_mc_lund.putByte("parent",    i_par, static_cast<int8_t>(par.mother1));
          b_mc_lund.putByte("daughter",  i_par, static_cast<int8_t>(par.daughter1));
          b_mc_lund.putFloat("px",       i_par, static_cast<float>(par.px));
          b_mc_lund.putFloat("py",       i_par, static_cast<float>(par.py));
          b_mc_lund.putFloat("pz",       i_par, static_cast<float>(par.pz));
          b_mc_lund.putFloat("energy",   i_par, static_cast<float>(par.energy));
          b_mc_lund.putFloat("mass",     i_par, static_cast<float>(par.mass));
          b_mc_lund.putFloat("vx",       i_par, static_cast<float>(par.vx));
          b_mc_lund.putFloat("vy",       i_par, static_cast<float>(par.vy));
          b_mc_lund.putFloat("vz",       i_par, static_cast<float>(par.vz));
          i_par++;
        }
        // fill `MC::User`
        b_mc_user.setRows(head.user_values.size());
        int i_user = 0;
        for(auto const& val : head.user_values)
          b_mc_user.putFloat("userVar", i_user++, static_cast<float>(val));
        // write
        event.addStructure(b_mc_header);
        event.addStructure(b_mc_event);
        event.addStructure(b_mc_particle);
        event.addStructure(b_mc_lund);
        event.addStructure(b_mc_user);
        writer.addEvent(event);
        // reset
        event.reset();
        b_mc_header.reset();
        b_mc_event.reset();
        b_mc_particle.reset();
        b_mc_lund.reset();
        b_mc_user.reset();
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
