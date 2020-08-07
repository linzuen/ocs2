#pragma once

#include <ocs2_cart_pole_example/definitions.h>
#include <ocs2_core/Types.h>

namespace ocs2 {
namespace cartpole {

/**
 * @brief Convert ocs2 cartpole state to generalized coordinate and generalized velocity used by RAIsim
 * @param[in] state: the state to be converted
 * @return {q, dq} pair that represents the state
 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> stateToRaisimGenCoordGenVel(const vector_t& state, const vector_t&) {
  Eigen::VectorXd q(2), dq(2);
  q << state(1), state(0);
  dq << state(3), state(2);
  return {q, dq};
}

/**
 * @brief Convert RAIsim generalized coordinates and velocities to ocs2 cartpole state
 * @note This should be the inverse to stateToRaisimGenCoordGenVel
 * @param[in] q the generalized coordinate
 * @param[in] dq the generalized velocity
 * @return the corresponding ocs2 cart pole state
 */
vector_t raisimGenCoordGenVelToState(const Eigen::VectorXd& q, const Eigen::VectorXd& dq) {
  vector_t state(STATE_DIM);
  state << q(1), q(0), dq(1), dq(0);
  return state;
}

/**
 * @brief Convert ocs2 control input to RAIsim generalized force
 * @param[in] input: The control computed by the ocs2 controller
 * @param[in] state: The current state
 * @return The generalized forces to be applied to the system
 */
Eigen::VectorXd inputToRaisimGeneralizedForce(double, const vector_t& input, const vector_t&, const Eigen::VectorXd&,
                                              const Eigen::VectorXd&) {
  Eigen::VectorXd generalizedForce = Eigen::VectorXd::Zero(2);
  generalizedForce.tail(1) = input;
  return generalizedForce;
}

}  // namespace cartpole
}  // namespace ocs2
