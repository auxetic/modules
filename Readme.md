**仅供参考**

### 下载
```bash
git clone git@git.ustclug.org:liuxu/modules.git     # 需要配置ssh-key
or
git clone https://git.ustclug.org/liuxu/modules.git

cd modules
```
### 使用
从复制 templates 目录拷贝示例主程序
```bash
cp templates/calc_msd.f90 00-main.f90
```
编译
```bash
mkdir build
cd build
cmake ..
make
```
