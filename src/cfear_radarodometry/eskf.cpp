#include "cfear_radarodometry/odometrykeyframefuser.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <chrono>
#include <thread>
#include <Sophus/so3.hpp>
#include <Eigen/Geometry>
#include <sophus/se3.hpp>


void eskf::eskf(pcl::Registration<PointT, PointT>::Ptr& registration, const ros::Time& stamp, const Eigen::Vector3f& pos, const Eigen::Quaternionf& quat)
{
  //单位阵初始化，随后给过程噪声。
  process_noise = Eigen::MatrixXf::Identity(18,18);
  process_noise.middleRows(0, 3) *= 1.0;
  process_noise.middleRows(3, 3) *= 1.0;
  process_noise.middleRows(6, 4) *= 0.5;
  process_noise.middleRows(10, 3) *= 1e-6;
  process_noise.middleRows(13, 3) *= 1e-6;

  //测量噪声，单位阵
  Eigen::MatrixXf measurement_noise = Eigen::MatrixXf::Identity(7, 7);
  measurement_noise.middleRows(0, 3) *= 0.01;
  measurement_noise.middleRows(3, 4) *= 0.001;

  //初始化协方差
  Eigen::MatrixXf cov = Eigen::MatrixXf::Identity(16, 16) * 0.01;

}

void eskf::Predict()(const ros::Time& stamp, const Eigen::Vector3f& acc, const Eigen::Vector3f& gyro) 
{

}

//状态方程f
VectorXt f(const VectorXt& state, const VectorXt& control) const {
    VectorXt next_state(16);

    Vector3t pt = state.middleRows(0, 3);  //位置
    Vector3t vt = state.middleRows(3, 3);  //速度
    Quaterniont qt(state[6], state[7], state[8], state[9]);
    qt.normalize(); // 归一化四元数

    Vector3t acc_bias = state.middleRows(10, 3);  //加速度偏差
    Vector3t gyro_bias = state.middleRows(13, 3); //陀螺仪偏差

    Vector3t raw_acc = control.middleRows(0, 3);  //加速度控制
    Vector3t raw_gyro = control.middleRows(3, 3);  //陀螺仪控制
	
	//下一时刻状态
    // position 首先更新位置
    next_state.middleRows(0, 3) = pt + vt * dt;

    // velocity。 更新速度，实际上并没有利用加速度矫正速度，原因是认为加速度噪声较大，对最终的精度并没有贡献。
    Vector3t g(0.0f, 0.0f, -9.80665f);
    Vector3t acc_ = raw_acc - acc_bias;
    Vector3t acc = qt * acc_;
    next_state.middleRows(3, 3) = vt; // + (acc - g) * dt;		// acceleration didn't contribute to accuracy due to large noise

    // orientation。首先完成了陀螺仪的增量计算并归一化（直接转化为四元数形式），将其转换为下一时刻的四元数。
    Vector3t gyro = raw_gyro - gyro_bias;
    Quaterniont dq(1, gyro[0] * dt / 2, gyro[1] * dt / 2, gyro[2] * dt / 2);
    dq.normalize();
    Quaterniont qt_ = (qt * dq).normalized();
    next_state.middleRows(6, 4) << qt_.w(), qt_.x(), qt_.y(), qt_.z();
	//将当前控制量传入下一时刻的状态向量。认为加速度和角速度上偏差不变
    next_state.middleRows(10, 3) = state.middleRows(10, 3);		// constant bias on acceleration
    next_state.middleRows(13, 3) = state.middleRows(13, 3);		// constant bias on angular velocity

    return next_state;
  }

  // 观测方程
  VectorXt h(const VectorXt& state) const {
    VectorXt observation(7);
    observation.middleRows(0, 3) = state.middleRows(0, 3);
    observation.middleRows(3, 4) = state.middleRows(6, 4).normalized();

    return observation;
  }



class eskf
{

public:
  eskf(pcl::Registration<PointT, PointT>::Ptr& registration, const ros::Time& stamp, const Eigen::Vector3d& pos, const Eigen::Quaternionf& quat)
  ~eskf();


private:
  Eigen::Vector2d p_ = Eigen::Vector2d::Zero();  // 初始位置为零向量
  Eigen::Vector2d v_ = Eigen::Vector2d::Zero();  // 初始速度为零向量


  Eigen::Matrix2d R_ = Eigen::Matrix2d::Identity();  // 单位矩阵，表示初始方向为单位阵
  Sophus::SO2d R_SO2 ;



  double bg_ = 0.0;  // 初始陀螺仪偏差为零
  Eigen::Vector2d ba_ = Eigen::Vector2d::Zero();  // 初始加速度计偏差为零向量
//   Eigen::Vector3d g_{0, 0, -9.8};  // 初始重力向量，朝向地球中心，大小为9.8m/s²
  Eigen::Matrix<double, 9, 1> dx_ = Eigen::Matrix<double, 9, 1>::Zero();  // 初始误差状态为零向量

  Eigen::Matrix<double, 9, 9> cov_ = Eigen::Matrix<double, 9, 9>::Identity();  // 协方差阵初始化为单位阵
  Eigen::Matrix<double, 9, 9> Q_ = Eigen::Matrix<double, 9, 9>::Zero();  // 噪声阵初始化为零矩阵
};

//这里的
bool eskf::Predict(const Eigen::Vector2d& gyro, const double acce, double timestamp) {
    assert(timestamp >= current_time_);

    double dt = timestamp - current_time_;
    if (dt > (5 * options_.imu_dt_) || dt < 0) {
        // 时间间隔不对，可能是第一个IMU数据，没有历史信息
        LOG(INFO) << "skip this imu because dt_ = " << dt;
        current_time_ = timestamp;
        return false;
    }

    //名义状态 p77 3.41
    // nominal state 递推
    Eigen::Matrix<double, 9, 1> new_p = p_ + v_ * dt + 0.5 * (R_ * (acce - ba_)) * dt * dt ;
    Eigen::Matrix<double, 9, 1> new_v = v_ + R_ * (acce - ba_) * dt ;
    Sophus::SO2d new_R = R_ * Sophus::SO2d::exp((gyro - bg_) * dt);

    R_ = new_R;
    v_ = new_v;
    p_ = new_p;
    // 其余状态维度不变

    // error state 递推
    // 计算运动过程雅可比矩阵 F，见(3.47)
    // F实际上是稀疏矩阵，也可以不用矩阵形式进行相乘而是写成散装形式，这里为了教学方便，使用矩阵形式
    Eigen::Matrix<double, 9, 9> F = Eigen::Matrix<double, 9, 9>::Identity();                                                 // 主对角线
    F.block<2, 2>(0, 2) = Eigen::Matrix<double, 2, 2>::Identity() * dt;                         // p 对 v
    F.block<2, 2>(2, 4) = -R_.matrix() * Sophus::SO2d::hat(acce - ba_) * dt;  // v对theta
    F.block<2, 2>(2, 7) = -R_.matrix() * dt;                             // v 对 ba
    // F.block<3, 3>(2, 9) = Eigen::Matrix<double, 3, 3>::Identity() * dt;                        // v 对 g //但对于2d，直接可以把g给消除掉
    F.block<2, 2>(4, 4) = Sophus::SO2d::exp(-(gyro - bg_) * dt).matrix();     // theta 对 theta
    F.block<2, 1>(4, 6) = dt;                        // theta 对 bg

    // mean and cov prediction
    dx_ = F * dx_;  // 这行其实没必要算，dx_在重置之后应该为零，因此这步可以跳过，但F需要参与Cov部分计算，所以保留
    cov_ = F * cov_.eval() * F.transpose() + Q_;
    current_time_ = timestamp;
    return true;
}

//得到更新后的 名义状态变量x 和 协方差矩阵P ，然后对P阵进行投影 ，并把误差项 dx 置零
//我看用的是odometrykeyframefuser用的是affine3d
bool eskf::ObserveSE2(const Sophus::SE3d& pose, double trans_noise, double ang_noise) {
    /// 既有旋转，也有平移
    /// 观测状态变量中的p, R，H为6x18，其余为零
    Eigen::Matrix<double, 4, 9> H = Eigen::Matrix<double, 4, 9>::Zero();
    H.block<2, 2>(0, 0) = Eigen::Matrix<double, 2, 2>::Identity();  // P部分
    H.block<2, 2>(2, 4) = Eigen::Matrix<double, 2, 2>::Identity();  // R部分（3.66)

    // 卡尔曼增益和更新过程
    Eigen::Matrix<double, 4, 1> noise_vec;
    noise_vec << trans_noise, trans_noise,  ang_noise, ang_noise;

    Eigen::Matrix<double, 4, 4> V = noise_vec.asDiagonal();
    Eigen::Matrix<double, 9, 4> K = cov_ * H.transpose() * (H * cov_ * H.transpose() + V).inverse();

    // 更新x和cov
    Eigen::Matrix<double, 4, 1> innov = Eigen::Matrix<double, 4, 1>::Zero();
    innov.head<2>() = (pose.translation() - p_);          // 平移部分
    innov.tail<2>() = (R_.inverse() * pose.so2()).log();  // 旋转部分(3.67)

    //innov是残差
    dx_ = K * innov;
    cov_ = (Eigen::Matrix<double, 9, 9>::Identity() - K * H) * cov_;

    UpdateAndReset();
    return true;
}

void eskf::UpdateAndReset() 
{
    p_ += dx_.block<2, 1>(0, 0);
    v_ += dx_.block<2, 1>(3, 0);
    R_ = R_ * Sophus::SO2d::exp(dx_.block<2, 1>(4, 0));

    if (options_.update_bias_gyro_) {
        bg_ += dx_.block<1, 1>(6, 0);
    }

    if (options_.update_bias_acce_) {
        ba_ += dx_.block<2, 1>(7, 0);
    }

    // g_ += dx_.block<3, 1>(15, 0);

    ProjectCov();
    dx_.setZero();
}

void eskf::ProjectCov() 
{
    Eigen::Matrix<double, 9, 9> J = Eigen::Matrix<double, 9, 9>::Identity();
    //dx 的R部分
    J.block<2, 2>(4, 4) = Eigen::Matrix<double, 2, 2>::Identity() - 0.5 * Sophus::SO2d::hat(dx_.block<2, 1>(4, 0));
    cov_ = J * cov_ * J.transpose();
}

void LooselyLIO::Predict() {
    imu_states_.clear();
    imu_states_.emplace_back(eskf_.GetNominalState()); //current_time_, R_, p_, v_, bg_, ba_

    /// 对IMU状态进行预测
    for (auto &imu : measures_.imu_) {
        eskf_.Predict(*imu);
        imu_states_.emplace_back(eskf_.GetNominalState());
    }
}


void run()
{
    LooselyLIO::Predict();
    /// 从EKF中获取预测pose，放入LO，获取LO位姿，最后合入EKF
    SE3 pose_predict = eskf_.GetNominalSE3();
    inc_ndt_lo_->AddCloud(current_scan_filter, pose_predict, true);
    pose_of_lo_ = pose_predict;
    eskf_.ObserveSE3(pose_of_lo_, 1e-2, 1e-2);
}


// 将Sophus::SO2转换为Eigen::Matrix2d
Eigen::Matrix2d SO2ToMatrix(const Sophus::SO2d& so2) {
    return so2.matrix();
}

// 将Eigen::Matrix2d转换为Sophus::SO2d
Sophus::SO2d MatrixToSO2(const Eigen::Matrix2d& matrix) {
    return Sophus::SO2d(matrix);
}


// 将Sophus::SO3转换为Eigen::Matrix3d
Eigen::Matrix3d SO3ToMatrix(const Sophus::SO3d& so3) {
    return so3.matrix();
}

// 将Eigen::Matrix3d转换为Sophus::SO3
Sophus::SO3d MatrixToSO3(const Eigen::Matrix3d& matrix) {
    return Sophus::SO3d(matrix);
}

// 将 Eigen::Affine3d 转换为 Sophus::SE3d
Sophus::SE3d affineToSE3(const Eigen::Affine3d& affine) {
    return Sophus::SE3d(affine.matrix());
}

// 将 Sophus::SE3d 转换为 Eigen::Affine3d
Eigen::Affine3d SE3ToAffine(const Sophus::SE3d& se3) {
    return Eigen::Affine3d(se3.matrix());
}

