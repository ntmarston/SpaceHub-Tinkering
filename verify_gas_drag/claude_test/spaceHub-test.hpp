// Modified spaceHub.hpp for isolation testing
// Includes disk-model-test.hpp (local) instead of disk-model.hpp (original)
// All other includes point to ../../SpaceHub/src/

#pragma once

#ifdef MPFR_VERSION_MAJOR
#include "../../SpaceHub/src/mpfr.hpp"
#endif

#include "../../SpaceHub/src/args-callback/callbacks.hpp"
#include "../../SpaceHub/src/args-callback/collision.hpp"
#include "../../SpaceHub/src/integrator/Gauss-Radau.hpp"
#include "../../SpaceHub/src/integrator/symplectic/symplectic-integrator.hpp"
#include "../../SpaceHub/src/interaction/alpha-disk.hpp"
#include "../../SpaceHub/src/interaction/drag-forces.hpp"
#include "disk-model-test.hpp"  // LOCAL MODIFIED VERSION
#include "../../SpaceHub/src/interaction/magneto-disk.hpp"
#include "../../SpaceHub/src/interaction/newtonian.hpp"
#include "../../SpaceHub/src/interaction/post-newtonian.hpp"
#include "../../SpaceHub/src/interaction/tidal.hpp"
#include "../../SpaceHub/src/kahan-number.hpp"
#include "../../SpaceHub/src/macros.hpp"
#include "../../SpaceHub/src/multi-thread/multi-thread.hpp"
#include "../../SpaceHub/src/ode-iterator/Bulirsch-Stoer.hpp"
#include "../../SpaceHub/src/ode-iterator/IAS15.hpp"
#include "../../SpaceHub/src/ode-iterator/const-iterator.hpp"
#include "../../SpaceHub/src/ode-iterator/error-checker/RMS.hpp"
#include "../../SpaceHub/src/ode-iterator/error-checker/max-ratio-error.hpp"
#include "../../SpaceHub/src/ode-iterator/error-checker/worst-offender.hpp"
#include "../../SpaceHub/src/ode-iterator/sequent-iterator.hpp"
#include "../../SpaceHub/src/ode-iterator/step-controller/PID-controller.hpp"
#include "../../SpaceHub/src/ode-iterator/step-controller/const-controller.hpp"
#include "../../SpaceHub/src/orbits/orbits.hpp"
#include "../../SpaceHub/src/orbits/particle-manip.hpp"
#include "../../SpaceHub/src/particle-system/archain.hpp"
#include "../../SpaceHub/src/particle-system/base-system.hpp"
#include "../../SpaceHub/src/particle-system/chain-system.hpp"
#include "../../SpaceHub/src/particle-system/regu-system.hpp"
#include "../../SpaceHub/src/particles/drag-particles.hpp"
#include "../../SpaceHub/src/particles/finite-size.hpp"
#include "../../SpaceHub/src/particles/point-particles.hpp"
#include "../../SpaceHub/src/particles/tide-particles.hpp"
#include "../../SpaceHub/src/scattering/cross-section.hpp"
#include "../../SpaceHub/src/scattering/hierarchical.hpp"
#include "../../SpaceHub/src/simulator.hpp"
#include "../../SpaceHub/src/stellar/stellar.hpp"
#include "../../SpaceHub/src/tools/auto-name.hpp"
#include "../../SpaceHub/src/tools/config-reader.hpp"
#include "../../SpaceHub/src/tools/timer.hpp"
#include "../../SpaceHub/src/type-class.hpp"

namespace hub {

#define USING_NAMESPACE_SPACEHUB_ALL \
    using namespace hub;             \
    using namespace hub::calc;       \
    using namespace hub::tools;      \
    using namespace hub::ode;        \
    using namespace hub::integrator; \
    using namespace hub::orbit;      \
    using namespace hub::unit;       \
    using namespace hub::particles;  \
    using namespace hub::random;     \
    using namespace hub::callback;   \
    using namespace hub::system;     \
    using namespace hub::force

    using DefaultTypes = Types<double, Vec3>;

    using DefaultForce = force::Interactions<hub::force::NewtonianGrav>;

    template <typename T>
    using DefaultParticles = particles::PointParticles<T>;
    namespace methods {
        namespace details {
            using namespace ode;
            using namespace integrator;
            using normal_type = Types<double, Vec3>;
            using extended_type = Types<long double, Vec3>;
            using precise_type = Types<double_k, Vec3>;
            using extended_precise_type = Types<long_double_k, Vec3>;
#ifdef MPFR_VERSION_MAJOR
            using any_bits_type = Types<mpfr::mpreal, Vec3>;
            using precise_any_bits_type = Types<mpreal_k, Vec3>;
#endif
            using rms_err = ode::RMS<normal_type>;
            using worst_offender_err = ode::WorstOffender<normal_type>;
            using adaptive_step_ctrl = PIDController<normal_type>;
            using const_step_ctrl = ConstStepController<normal_type>;

            using rms_err_ext = ode::RMS<extended_type>;
            using worst_offender_err_ext = ode::WorstOffender<extended_type>;
            using adaptive_step_ctrl_ext = PIDController<extended_type>;
            using const_step_ctrl_ext = ConstStepController<extended_type>;

            using const_sym2 = ConstOdeIterator<Symplectic2nd<normal_type>>;
            using const_sym4 = ConstOdeIterator<Symplectic4th<normal_type>>;
            using const_sym6 = ConstOdeIterator<Symplectic6th<normal_type>>;
            using const_sym8 = ConstOdeIterator<Symplectic8th<normal_type>>;
            using const_sym10 = ConstOdeIterator<Symplectic10th<normal_type>>;
            using const_Radau = ConstOdeIterator<GaussRadau<normal_type>>;

            using const_sym2_ext = ConstOdeIterator<Symplectic2nd<extended_type>>;
            using const_sym4_ext = ConstOdeIterator<Symplectic4th<extended_type>>;
            using const_sym6_ext = ConstOdeIterator<Symplectic6th<extended_type>>;
            using const_sym8_ext = ConstOdeIterator<Symplectic8th<extended_type>>;
            using const_sym10_ext = ConstOdeIterator<Symplectic10th<extended_type>>;
            using const_Radau_ext = ConstOdeIterator<GaussRadau<extended_type>>;

            using const_sym2_plus = ConstOdeIterator<Symplectic2nd<precise_type>>;
            using const_sym4_plus = ConstOdeIterator<Symplectic4th<precise_type>>;
            using const_sym6_plus = ConstOdeIterator<Symplectic6th<precise_type>>;
            using const_sym8_plus = ConstOdeIterator<Symplectic8th<precise_type>>;
            using const_sym10_plus = ConstOdeIterator<Symplectic10th<precise_type>>;
            using const_Radau_plus = ConstOdeIterator<GaussRadau<precise_type>>;

            using const_sym2_extplus = ConstOdeIterator<Symplectic2nd<extended_precise_type>>;
            using const_sym4_extplus = ConstOdeIterator<Symplectic4th<extended_precise_type>>;
            using const_sym6_extplus = ConstOdeIterator<Symplectic6th<extended_precise_type>>;
            using const_sym8_extplus = ConstOdeIterator<Symplectic8th<extended_precise_type>>;
            using const_sym10_extplus = ConstOdeIterator<Symplectic10th<extended_precise_type>>;
            using const_Radau_extplus = ConstOdeIterator<GaussRadau<extended_precise_type>>;

            using BS = BulirschStoer<LeapFrogDKD<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym2 = SequentOdeIterator<Symplectic2nd<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym4 = SequentOdeIterator<Symplectic4th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym6 = SequentOdeIterator<Symplectic6th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym8 = SequentOdeIterator<Symplectic8th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym10 = SequentOdeIterator<Symplectic10th<normal_type>, worst_offender_err, adaptive_step_ctrl>;
            using Radau = IAS15<GaussRadau<normal_type>, MaxRatioError<normal_type>, adaptive_step_ctrl>;

            using BS_ext = BulirschStoer<LeapFrogDKD<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym2_ext =
                SequentOdeIterator<Symplectic2nd<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym4_ext =
                SequentOdeIterator<Symplectic4th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym6_ext =
                SequentOdeIterator<Symplectic6th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym8_ext =
                SequentOdeIterator<Symplectic8th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym10_ext =
                SequentOdeIterator<Symplectic10th<extended_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using Radau_ext = IAS15<GaussRadau<extended_type>, MaxRatioError<extended_type>, adaptive_step_ctrl_ext>;

            using BS_plus = BulirschStoer<LeapFrogDKD<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym2_plus = SequentOdeIterator<Symplectic2nd<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym4_plus = SequentOdeIterator<Symplectic4th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym6_plus = SequentOdeIterator<Symplectic6th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym8_plus = SequentOdeIterator<Symplectic8th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using sym10_plus = SequentOdeIterator<Symplectic10th<precise_type>, worst_offender_err, adaptive_step_ctrl>;
            using Radau_plus = IAS15<GaussRadau<precise_type>, MaxRatioError<normal_type>, adaptive_step_ctrl>;

            using BS_extplus =
                BulirschStoer<LeapFrogDKD<extended_precise_type>, worst_offender_err_ext, adaptive_step_ctrl_ext>;
            using sym2_extplus = SequentOdeIterator<Symplectic2nd<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym4_extplus = SequentOdeIterator<Symplectic4th<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym6_extplus = SequentOdeIterator<Symplectic6th<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym8_extplus = SequentOdeIterator<Symplectic8th<extended_precise_type>, worst_offender_err_ext,
                                                    adaptive_step_ctrl_ext>;
            using sym10_extplus = SequentOdeIterator<Symplectic10th<extended_precise_type>, worst_offender_err_ext,
                                                     adaptive_step_ctrl_ext>;
            using Radau_extplus =
                IAS15<GaussRadau<extended_precise_type>, MaxRatioError<extended_type>, adaptive_step_ctrl_ext>;
#ifdef MPFR_VERSION_MAJOR
            using ABits = BulirschStoer<LeapFrogDKD<any_bits_type>, ode::WorstOffender<any_bits_type>,
                                        PIDController<any_bits_type>, 32>;
            using ABits_plus =
                BulirschStoer<LeapFrogDKD<precise_any_bits_type>, ode::WorstOffender<precise_any_bits_type>,
                              PIDController<precise_any_bits_type>, 32>;
#endif
        };  // namespace details

#define DEFINE_ADAPTIVE_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                         \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME = Simulator<system::SYSTEM<particle<details::normal_type>, interactions>, details::ITER>;               \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Plus =                                                                                                \
        Simulator<system::SYSTEM<particle<details::precise_type>, interactions>, details::ITER##_plus>;                \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Ext = Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::ITER##_ext>; \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_ExtPlus =                                                                                             \
        Simulator<system::SYSTEM<particle<details::extended_precise_type>, interactions>, details::ITER##_extplus>;

#define DEFINE_CONST_STEP_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME =                                                                                          \
        Simulator<particle_system::SYSTEM<particle<details::normal_type>, interactions>, details::const_##ITER>;  \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_Plus = Simulator<particle_system::SYSTEM<particle<details::precise_type>, interactions>, \
                                          details::const_##ITER##_plus>;                                          \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_Ext = Simulator<particle_system::SYSTEM<particle<details::extended_type>, interactions>, \
                                         details::const_##ITER##_ext>;                                            \
                                                                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>     \
    using Const_##NAME##_ExtPlus =                                                                                \
        Simulator<particle_system::SYSTEM<particle<details::extended_precise_type>, interactions>,                \
                  details::const_##ITER##_extplus>;

#define DEFINE_INTEGRATION_METHOD(NAME, SYSTEM, ITER)                                                                  \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME = Simulator<system::SYSTEM<particle<details::normal_type>, interactions>, details::ITER>;               \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Plus =                                                                                                \
        Simulator<system::SYSTEM<particle<details::precise_type>, interactions>, details::ITER##_plus>;                \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME =                                                                                               \
        Simulator<system::SYSTEM<particle<details::normal_type>, interactions>, details::const_##ITER>;                \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME##_Plus =                                                                                        \
        Simulator<system::SYSTEM<particle<details::precise_type>, interactions>, details::const_##ITER##_plus>;        \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_Ext = Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::ITER##_ext>; \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using NAME##_ExtPlus =                                                                                             \
        Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::ITER##_extplus>;            \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME##_Ext =                                                                                         \
        Simulator<system::SYSTEM<particle<details::extended_type>, interactions>, details::const_##ITER##_ext>;        \
                                                                                                                       \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>          \
    using Const_##NAME##_ExtPlus = Simulator<system::SYSTEM<particle<details::extended_precise_type>, interactions>,   \
                                             details::const_##ITER##_extplus>;

#define DEFINE_ADAPTIVE_ARBITRARY_BIT_METHOD(NAME, SYSTEM, ITER)                                              \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles> \
    using NAME = Simulator<system::SYSTEM<particle<details::any_bits_type>, interactions>, details::ITER>;    \
                                                                                                              \
    template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles> \
    using NAME##_Plus =                                                                                       \
        Simulator<system::SYSTEM<particle<details::precise_any_bits_type>, interactions>, details::ITER##_plus>;

        DEFINE_ADAPTIVE_INTEGRATION_METHOD(BS, SimpleSystem, BS)
        DEFINE_ADAPTIVE_INTEGRATION_METHOD(AR_BS, RegularizedSystem, BS)
        DEFINE_ADAPTIVE_INTEGRATION_METHOD(Chain_BS, ChainSystem, BS)
        DEFINE_ADAPTIVE_INTEGRATION_METHOD(AR_Chain, ARchainSystem, BS)
#ifdef MPFR_VERSION_MAJOR
        DEFINE_ADAPTIVE_ARBITRARY_BIT_METHOD(ABITS, SimpleSystem, ABits)
        DEFINE_ADAPTIVE_ARBITRARY_BIT_METHOD(AR_ABITS, RegularizedSystem, ABits)
#endif
        DEFINE_INTEGRATION_METHOD(Sym2, SimpleSystem, sym2)
        DEFINE_INTEGRATION_METHOD(AR_Sym2, RegularizedSystem, sym2)
        DEFINE_INTEGRATION_METHOD(Chain_Sym2, ChainSystem, sym2)
        DEFINE_INTEGRATION_METHOD(AR_Sym2_Chain, ARchainSystem, sym2)
        DEFINE_INTEGRATION_METHOD(Sym4, SimpleSystem, sym4)
        DEFINE_INTEGRATION_METHOD(AR_Sym4, RegularizedSystem, sym4)
        DEFINE_INTEGRATION_METHOD(Chain_Sym4, ChainSystem, sym4)
        DEFINE_INTEGRATION_METHOD(AR_Sym4_Chain, ARchainSystem, sym4)
        DEFINE_INTEGRATION_METHOD(Sym6, SimpleSystem, sym6)
        DEFINE_INTEGRATION_METHOD(AR_Sym6, RegularizedSystem, sym6)
        DEFINE_INTEGRATION_METHOD(Chain_Sym6, ChainSystem, sym6)
        DEFINE_INTEGRATION_METHOD(AR_Sym6_Chain, ARchainSystem, sym6)
        DEFINE_INTEGRATION_METHOD(Sym8, SimpleSystem, sym8)
        DEFINE_INTEGRATION_METHOD(AR_Sym8, RegularizedSystem, sym8)
        DEFINE_INTEGRATION_METHOD(Chain_Sym8, ChainSystem, sym8)
        DEFINE_INTEGRATION_METHOD(AR_Sym8_Chain, ARchainSystem, sym8)
        DEFINE_INTEGRATION_METHOD(Sym10, SimpleSystem, sym10)
        DEFINE_INTEGRATION_METHOD(AR_Sym10, RegularizedSystem, sym10)
        DEFINE_INTEGRATION_METHOD(Chain_Sym10, ChainSystem, sym10)
        DEFINE_INTEGRATION_METHOD(AR_Sym10_Chain, ARchainSystem, sym10)
        DEFINE_INTEGRATION_METHOD(Radau, SimpleSystem, Radau)
        DEFINE_INTEGRATION_METHOD(AR_Radau, RegularizedSystem, Radau)
        DEFINE_INTEGRATION_METHOD(Chain_Radau, ChainSystem, Radau)
        DEFINE_INTEGRATION_METHOD(AR_Radau_Chain, ARchainSystem, Radau)

        template <typename interactions = DefaultForce, template <typename> typename particle = DefaultParticles>
        using DefaultMethod = methods::AR_Chain_Plus<interactions, particle>;
    }  // namespace methods

    template <typename T>
    inline constexpr void set_mpreal_bits_from_rtol(T rtol) {
#ifdef MPFR_VERSION_MAJOR
        mpfr::mpreal::set_default_prec(size_t(mpfr::fabs(mpfr::LOG10(rtol))) * 4 + 32);
#endif
    }
}  // namespace hub
