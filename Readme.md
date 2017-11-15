**说明**
------

### 提示
* 需要使用此代码进行工作的话, 建议fork此库，并在fork后的代码库中，新建work分支。
```bash
点击上方Fork按钮，转到自己fork源的页面
git clone git@git.ustclug.org:"""your user name in git.ustclug"""/modules.git
cd modules
git checkout -b work
```
* fork后的代码和上游同步。
```bash
git remote add upstream https://git.ustclug.org/liuxu/modules.git
git checkout master
git pull upstream master
git push 
```
* 更新work分支，master与上游同步后，使用变基(rebase)同步work分支
```bash
git checkout work
git rebase -i master
解决冲突，合并
```
* 欢迎merge request

### 下载
```bash
git clone git@git.ustclug.org:liuxu/modules.git     # 需要配置ssh-key
or
git clone https://git.ustclug.org/liuxu/modules.git
cd modules
```
### 使用
从 templates 目录拷贝示例主程序
```bash
cp templates/calc_msd.f90 main.f90
```
编译
```bash
mkdir build && cd build
export FC=ifort     # 可选，gfortran 版本过低（4.4）时，编译会报错。
cmake ..
make                # ubuntu可能链接lapack失败，原因未知，请使用RedHat系linux系统
./main
```


### 模块说明
-----
#### syst
1. 体系维度
2. 大小球比例
...


#### config
1. 定义tpcon结构，包含体系粒子数、位形、半径、速度等等，盒子边长、剪切等等，系统能量、压强等等。
2. 定义用于初始化系统的子过程：init_system
3. 用于产生位形的：gen_rand_config, gen_lattice_triangle
4. 用于计算两粒子距离的子过程：calc_dra
5. 用于裁剪位形的子过程：trim_config


#### list
1. 定义tplist、tplistone结构，包含粒子邻居编号，ra_ij除以盒子边长后的取整信息
2. 定义子过程分别对应列表空间分配、列表、列表检查

* [ ] 检查自由粒子（rattler）
* [ ] 配位数


#### fire
1. 子过程init_confire，复制欲minimize的体统
2. 子过程check_system_force，检查体系残余力
3. 子过程mini_fire_cv，恒定体积fire
4. 子过程mini_fire_cp，恒定压强和剪切力fire，可分别或控制xy边界


#### md
1. 分子动力学相关
2. 子过程md_nvt，实现NVT MD模拟的一步积分

* [ ] npt


#### force
1. 子过程calc_foce计算作用力，需要列表信息
2. 子过程calc_force_withoutlist，无需列表

* [ ] Lenard-Jones


#### network
1. 定义tpspring，tpnetwork结构，记录弹簧网络体系
2. 子过程make_network，从粒子位形和接触信息构造弹簧网络
3. 子过程calc_len，计算弹簧网络某一bond的长度


#### dynamic
1. 动力学量的计算
2. msd, fkt, 速度关联


#### ToDo struc
1. 静态结构计算

* [ ] gr sk
* [ ] psi4 psi6


#### matrix
1. 动力学矩阵，振动模相关
2. 子过程make_dymatrix，计算动力学矩阵
3. 子过程solve_matrix，调用Lapack解矩阵本征值，可选择只解前rangevar部分
4. 子过程calc_pw，计算参与度

* [ ] 考虑盒子可压缩膨胀和剪切的情形
* [ ] Psi_th


#### ToDo Voronois Cell
