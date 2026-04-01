#pragma once

#include "Lund.h"

#ifdef CLAS_STRINGSPINNER_HIPO
#include <hipo4/writer.h>

namespace css {

  class Hipo {

    private:
      hipo::writer writer;
      hipo::event event;
      hipo::schema s_mc_particle{"MC::Particle", 40, 2};
      hipo::bank b_mc_particle;

    public: // ==================================================================================

      Hipo(std::string hipo_file_name)
      {
        s_mc_particle.parse("pid/I,px/F,py/F,pz/F,vx/F,vy/F,vz/F,vt/F");
        writer.getDictionary().addSchema(s_mc_particle);
        b_mc_particle = hipo::bank(s_mc_particle);
        writer.open(hipo_file_name.c_str());
      }

      // ==================================================================================

      void SetNumParticles(int const num)
      {
        b_mc_particle.setRows(num);
      }

      // ==================================================================================

      void FillParticle(LundParticle const& lund, int const& row)
      {
        b_mc_particle.putInt("pid", row, lund.pdg);
        b_mc_particle.putFloat("px", row, lund.px);
        b_mc_particle.putFloat("py", row, lund.py);
        b_mc_particle.putFloat("pz", row, lund.pz);
        b_mc_particle.putFloat("vx", row, lund.vx);
        b_mc_particle.putFloat("vy", row, lund.vy);
        b_mc_particle.putFloat("vz", row, lund.vz);
        b_mc_particle.putFloat("vt", row, 0);
      }

      // ==================================================================================

      void Write()
      {
        event.reset();
        event.addStructure(b_mc_particle);
        writer.addEvent(event);
        b_mc_particle.reset();
      }

      // ==================================================================================

      void Close()
      {
        writer.close();
      }

  };
}
#endif
