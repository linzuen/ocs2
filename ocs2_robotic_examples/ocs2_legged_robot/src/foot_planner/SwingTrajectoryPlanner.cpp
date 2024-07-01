/******************************************************************************
Copyright (c) 2021, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

 * Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "ocs2_legged_robot/foot_planner/SwingTrajectoryPlanner.h"

#include <ocs2_core/misc/Lookup.h>

#include "ocs2_legged_robot/gait/MotionPhaseDefinition.h"

namespace ocs2 {
namespace legged_robot {

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
SwingTrajectoryPlanner::SwingTrajectoryPlanner(Config config, size_t numFeet) : config_(std::move(config)), numFeet_(numFeet) {}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t SwingTrajectoryPlanner::getZvelocityConstraint(size_t leg, scalar_t time) const {
  const auto index = lookup::findIndexInTimeArray(feetHeightTrajectoriesEvents_[leg], time);
  return feetHeightTrajectories_[leg][index].velocity(time);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t SwingTrajectoryPlanner::getZpositionConstraint(size_t leg, scalar_t time) const {
  const auto index = lookup::findIndexInTimeArray(feetHeightTrajectoriesEvents_[leg], time);
  return feetHeightTrajectories_[leg][index].position(time);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void SwingTrajectoryPlanner::update(const ModeSchedule& modeSchedule, scalar_t terrainHeight) {
  const scalar_array_t terrainHeightSequence(modeSchedule.modeSequence.size(), terrainHeight);
  feet_array_t<scalar_array_t> liftOffHeightSequence;
  liftOffHeightSequence.fill(terrainHeightSequence);
  feet_array_t<scalar_array_t> touchDownHeightSequence;
  touchDownHeightSequence.fill(terrainHeightSequence);
  update(modeSchedule, liftOffHeightSequence, touchDownHeightSequence);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void SwingTrajectoryPlanner::update(const ModeSchedule& modeSchedule, const feet_array_t<scalar_array_t>& liftOffHeightSequence,
                                    const feet_array_t<scalar_array_t>& touchDownHeightSequence) {
  const auto& modeSequence = modeSchedule.modeSequence;
  const auto& eventTimes = modeSchedule.eventTimes;

    // todo 返回4条腿 根据modeSchedule.modeSequence中的每个modeNumber转换为对应的contactFlag true or false
  const auto eesContactFlagStocks = extractContactFlags(modeSequence);

  feet_array_t<std::vector<int>> startTimesIndices;
  feet_array_t<std::vector<int>> finalTimesIndices;
    // * for 循环4条腿 ，返回每条腿的swing leg项 开始的index和结束的index，站立腿的index为0
  for (size_t leg = 0; leg < numFeet_; leg++) {
    std::tie(startTimesIndices[leg], finalTimesIndices[leg]) = updateFootSchedule(eesContactFlagStocks[leg]);
  }
    // 也是4条腿的循环 
  for (size_t j = 0; j < numFeet_; j++) {
    feetHeightTrajectories_[j].clear();
    feetHeightTrajectories_[j].reserve(modeSequence.size());
    for (int p = 0; p < modeSequence.size(); ++p) {
      // INFO 如果是swing leg
      if (!eesContactFlagStocks[j][p]) {  // for a swing leg
        const int swingStartIndex = startTimesIndices[j][p];
        const int swingFinalIndex = finalTimesIndices[j][p];
                // * 只是判断index是否为-1 或者大于等于modeSequence.size()，如果是则抛出异常
        checkThatIndicesAreValid(j, p, swingStartIndex, swingFinalIndex, modeSequence);

          // * 找到对应的swing leg的 startTime和finalTime，其实就是对应的modeSequence的时间，切换mode时间
        const scalar_t swingStartTime = eventTimes[swingStartIndex];
        const scalar_t swingFinalTime = eventTimes[swingFinalIndex];

        const scalar_t scaling = swingTrajectoryScaling(swingStartTime, swingFinalTime, config_.swingTimeScale);

        const CubicSpline::Node liftOff{swingStartTime, liftOffHeightSequence[j][p], scaling * config_.liftOffVelocity};
        const CubicSpline::Node touchDown{swingFinalTime, touchDownHeightSequence[j][p], scaling * config_.touchDownVelocity};
        const scalar_t midHeight = std::min(liftOffHeightSequence[j][p], touchDownHeightSequence[j][p]) + scaling * config_.swingHeight;
        feetHeightTrajectories_[j].emplace_back(liftOff, midHeight, touchDown);
      } else {  // for a stance leg
        // Note: setting the time here arbitrarily to 0.0 -> 1.0 makes the assert in CubicSpline fail
        const CubicSpline::Node liftOff{0.0, liftOffHeightSequence[j][p], 0.0};
        const CubicSpline::Node touchDown{1.0, liftOffHeightSequence[j][p], 0.0};
        feetHeightTrajectories_[j].emplace_back(liftOff, liftOffHeightSequence[j][p], touchDown);
      }
    }
    feetHeightTrajectoriesEvents_[j] = eventTimes;
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::pair<std::vector<int>, std::vector<int>> SwingTrajectoryPlanner::updateFootSchedule(const std::vector<bool>& contactFlagStock) {
    // info 就是modeSchedule.modeSequence中的size
  const size_t numPhases = contactFlagStock.size();

  //  todo startTimeIndexStock初始化为0，finalTimeIndexStock初始化为0
  std::vector<int> startTimeIndexStock(numPhases, 0);
  std::vector<int> finalTimeIndexStock(numPhases, 0);

  // find the startTime and finalTime indices for swing feet
  for (size_t i = 0; i < numPhases; i++) {
        // * contactFlagStock 为 false的时候表示是swing leg
    if (!contactFlagStock[i]) {
        // * 找到每个swing leg的startTimeIndexStock和finalTimeIndexStock，
      // info 比如 swing leg 的 startindex 为上一个 stance leg 的 index，finalindex 为下一个 stance leg 的 index
      std::tie(startTimeIndexStock[i], finalTimeIndexStock[i]) = findIndex(i, contactFlagStock);
    }
  }
  return {startTimeIndexStock, finalTimeIndexStock};
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
feet_array_t<std::vector<bool>> SwingTrajectoryPlanner::extractContactFlags(const std::vector<size_t>& phaseIDsStock) const {
    // info phaseIDsStock 是 modeSchedule.modeSequence;
  const size_t numPhases = phaseIDsStock.size();

  feet_array_t<std::vector<bool>> contactFlagStock;
    // todo 将4条腿的contactFlagStock初始化为4个空的vector<bool>，每个vector<bool>的长度为numPhases
  std::fill(contactFlagStock.begin(), contactFlagStock.end(), std::vector<bool>(numPhases));

  for (size_t i = 0; i < numPhases; i++) {
        // todo 将modeSchedule.modeSequence中的每个modeNumber转换为对应的contactFlag，返回 contact_flag_t
    const auto contactFlag = modeNumber2StanceLeg(phaseIDsStock[i]);
    for (size_t j = 0; j < numFeet_; j++) {
      contactFlagStock[j][i] = contactFlag[j];
    }
  }
    // todo 返回4条腿 根据modeSchedule.modeSequence中的每个modeNumber转换为对应的contactFlag true or false
  return contactFlagStock;
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
std::pair<int, int> SwingTrajectoryPlanner::findIndex(size_t index, const std::vector<bool>& contactFlagStock) {
  const size_t numPhases = contactFlagStock.size();

  // skip if it is a stance leg
    // * updateFootSchedule函数里面调用了findIndex
  // * 说明进入时候的判断失误了，在这里是不会进入的，估计是为了findIndex函数的通用性
  if (contactFlagStock[index]) {
    return {0, 0};
  }

  // find the starting time
  int startTimesIndex = -1;
  for (int ip = index - 1; ip >= 0; ip--) {
    if (contactFlagStock[ip]) {
      startTimesIndex = ip;
      break;
    }
  }

  // find the final time
  int finalTimesIndex = numPhases - 1;
  for (size_t ip = index + 1; ip < numPhases; ip++) {
    if (contactFlagStock[ip]) {
      finalTimesIndex = ip - 1;
      break;
    }
  }

  return {startTimesIndex, finalTimesIndex};
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
void SwingTrajectoryPlanner::checkThatIndicesAreValid(int leg, int index, int startIndex, int finalIndex,
                                                      const std::vector<size_t>& phaseIDsStock) {
  const size_t numSubsystems = phaseIDsStock.size();
  if (startIndex < 0) {
    std::cerr << "Subsystem: " << index << " out of " << numSubsystems - 1 << std::endl;
    for (size_t i = 0; i < numSubsystems; i++) {
      std::cerr << "[" << i << "]: " << phaseIDsStock[i] << ",  ";
    }
    std::cerr << std::endl;

    throw std::runtime_error("The time of take-off for the first swing of the EE with ID " + std::to_string(leg) + " is not defined.");
  }
  if (finalIndex >= numSubsystems - 1) {
    std::cerr << "Subsystem: " << index << " out of " << numSubsystems - 1 << std::endl;
    for (size_t i = 0; i < numSubsystems; i++) {
      std::cerr << "[" << i << "]: " << phaseIDsStock[i] << ",  ";
    }
    std::cerr << std::endl;

    throw std::runtime_error("The time of touch-down for the last swing of the EE with ID " + std::to_string(leg) + " is not defined.");
  }
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
scalar_t SwingTrajectoryPlanner::swingTrajectoryScaling(scalar_t startTime, scalar_t finalTime, scalar_t swingTimeScale) {
  return std::min(1.0, (finalTime - startTime) / swingTimeScale);
}

/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/
SwingTrajectoryPlanner::Config loadSwingTrajectorySettings(const std::string& fileName, const std::string& fieldName, bool verbose) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_info(fileName, pt);

  if (verbose) {
    std::cerr << "\n #### Swing Trajectory Config:";
    std::cerr << "\n #### =============================================================================\n";
  }

  SwingTrajectoryPlanner::Config config;
  const std::string prefix = fieldName + ".";

  loadData::loadPtreeValue(pt, config.liftOffVelocity, prefix + "liftOffVelocity", verbose);
  loadData::loadPtreeValue(pt, config.touchDownVelocity, prefix + "touchDownVelocity", verbose);
  loadData::loadPtreeValue(pt, config.swingHeight, prefix + "swingHeight", verbose);
  loadData::loadPtreeValue(pt, config.swingTimeScale, prefix + "swingTimeScale", verbose);

  if (verbose) {
    std::cerr << " #### =============================================================================" << std::endl;
  }

  return config;
}

}  // namespace legged_robot
}  // namespace ocs2
