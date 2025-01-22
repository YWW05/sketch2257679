// This is a remix of a beautiful sketch by @noel:
// https://openprocessing.org/sketch/872516
// My only changes were:
// 1) Drawing the particles in a circular (more cup-like?) enclosure.
// 2) Respawning them near the center of the circle, rather than at random x and y positions.
// 3) Gradually changing the color of the liquid and the particles.
// 4) Hiding the mouse pointer.
// 5) Adding a cup handle. 
// To be clear: I have no idea how the rest of the particle code works! :)
// - Ivan Rudnicki

// I added code to gradually change the mouseX and mouseY values. I then slowed
// the dissolve and color change to give time for the random movements to develop.
// Confession: I too have no idea how the rest of the particle code works! ;)
// - Richard Bourne

/**
 * Java implementation of the Navier-Stokes-Solver from
 * http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
 */

import KinectPV2.*;

KinectPV2 kinect;
NavierStokesSolver fluidSolver;
double visc, diff, vScale, velocityScale;
float limitVelocity;
float changeX;
float changeY;
int numParticles;
Particle[] particles;
int c = color(255);
color c1, c2, c3, surface;
boolean vectors = false;

// 调整深度阈值和其他参数
float minDepth = 500;    // 最小深度值
float maxDepth = 1000;   // 最大深度值
float depthScale = 2.0;  // 增加深度影响力度
float baseFlowSpeed = 0.5; // 基础流动速度
float flowNoiseScale = 0.002; // 噪声缩放因子
float flowNoiseTime = 0; // 噪声时间变量

// 添加前一帧的深度数据存储
int[] prevDepth;
float motionThreshold = 50; // 运动检测阈值

public void setup() {
  fullScreen();  // 使用全屏模式
  
  // 初始化 Kinect
  kinect = new KinectPV2(this);
  kinect.enableDepthImg(true);
  kinect.init();
  
  // 初始化前一帧深度数据
  prevDepth = new int[512 * 424];
  
  // 初始化颜色
  c1 = color(255, 0, 255);     // 亮粉色
  c2 = color(65, 105, 225);    // 皇家蓝
  c3 = color(0);               // 黑色
  
  // 调整流体参数
  fluidSolver = new NavierStokesSolver();
  numParticles = (int)pow(2, 13);  // 减少粒子数量
  particles = new Particle[numParticles];
  visc = 0.0001f;             // 降低粘性
  diff = 0.01f;               // 降低扩散率
  velocityScale = 30;         // 增加速度
  vectors = true;
  limitVelocity = 200;        // 增加最大速度限制
  initParticles();
  
  background(0);  // 初始黑色背景
  rectMode(CORNER);  // 确保rect从左上角开始绘制
}

private void initParticles() {
  for (int i = 0; i < numParticles - 1; i++) {
    particles[i] = new Particle();
    // 随机分布在整个屏幕范围内，而不是限制在圆形区域
    particles[i].x = random(width);
    particles[i].y = random(height);
  }
}

public void draw() {
  // 增加透明度使运动更明显
  fill(0, 35);  // 增加透明度
  rect(0, 0, width, height);
  
  // 获取深度数据并处理
  int[] rawDepth = kinect.getRawDepthData();
  processKinectData(rawDepth);
  
  // 更新流体
  double dt = 1 / frameRate;
  fluidSolver.tick(dt, visc, diff);
  vScale = velocityScale * 60. / frameRate;
  drawParticlesPixels();
}

void processKinectData(int[] rawDepth) {
  boolean motionDetected = false;
  int skip = 10;
  
  // Kinect深度图像分辨率
  int kw = 512;
  int kh = 424;
  
  for (int x = 0; x < kw; x += skip) {
    for (int y = 0; y < kh; y += skip) {
      int index = x + y * kw;
      if (index >= rawDepth.length) continue;
      
      float depth = rawDepth[index];
      float prevDepthValue = prevDepth[index];
      
      // 计算深度变化
      float depthChange = abs(depth - prevDepthValue);
      
      // 检测到显著运动
      if (depthChange > motionThreshold && depth > minDepth && depth < maxDepth) {
        motionDetected = true;
        float screenX = map(x, 0, kw, 0, width);
        float screenY = map(y, 0, kh, 0, height);
        
        // 计算运动方向和强度
        float dirX = (depth < prevDepthValue) ? 1 : -1;
        float dirY = (y > kh/2) ? 1 : -1;
        
        float force = pow(map(depthChange, 0, 500, 0, 1), 2) * depthScale;
        
        // 应用定向力
        float dx = dirX * force * 3;
        float dy = dirY * force * 3;
        
        int cellX = floor(screenX / (width / NavierStokesSolver.N));
        int cellY = floor(screenY / (height / NavierStokesSolver.N));
        fluidSolver.applyForce(cellX, cellY, dx, dy);
      }
    }
  }
  
  // 如果没有检测到运动，应用基础流动
  if (!motionDetected) {
    applyBaseFlow();
  }
  
  // 更新前一帧的深度数据
  arraycopy(rawDepth, 0, prevDepth, 0, rawDepth.length);
}

void applyBaseFlow() {
  // 使用Perlin噪声创建平滑的基础流动
  float noiseScale = flowNoiseScale;
  flowNoiseTime += 0.01;
  
  for (int i = 0; i < NavierStokesSolver.N; i++) {
    for (int j = 0; j < NavierStokesSolver.N; j++) {
      float angle = noise(i * noiseScale, j * noiseScale, flowNoiseTime) * TWO_PI * 2;
      float dx = cos(angle) * baseFlowSpeed;
      float dy = sin(angle) * baseFlowSpeed;
      fluidSolver.applyForce(i, j, dx, dy);
    }
  }
}

private void drawParticlesPixels() {
  int n = NavierStokesSolver.N;
  float cellHeight = height / n;
  float cellWidth = width / n;
  for (Particle p : particles) {
    if (p != null) {
      int cellX = floor(p.x / cellWidth);
      int cellY = floor(p.y / cellHeight);
      float dx = (float) fluidSolver.getDx(cellX, cellY);
      float dy = (float) fluidSolver.getDy(cellX, cellY);
      
      // 添加平滑的阻尼效果
      dx *= 0.99;
      dy *= 0.99;
      
      p.x += dx * vScale;
      p.y += dy * vScale;
			
      if (p.x < 0 || p.x > width || p.y < 0 || p.y > height) {
        p.x = random(width);
        p.y = random(height);
      }
      
      // 基于粒子位置计算颜色
      float colorMix = (sin(p.x * 0.01) + cos(p.y * 0.01)) * 0.5 + 0.5;
      color particleColor = lerpColor(c1, c2, colorMix);
      
      stroke(particleColor);
      strokeWeight(3);
      point(p.x, p.y);
    }
  }
}

public class Particle {
  public float x;
  public float y;
}

public class NavierStokesSolver {
  final static int N = 16;
  final static int SIZE = (N + 2) * (N + 2);
  double[] u = new double[SIZE];
  double[] v = new double[SIZE];
  double[] u_prev = new double[SIZE];
  double[] v_prev = new double[SIZE];
  public NavierStokesSolver() {
  }

  public double getDx(int x, int y) {
    return u[INDEX(x + 1, y + 1)];
  }

  public double getDy(int x, int y) {
    return v[INDEX(x + 1, y + 1)];
  }

  public void applyForce(int cellX, int cellY, double vx, double vy) {
    cellX += 1;
    cellY += 1;
    double dx = u[INDEX(cellX, cellY)];
    double dy = v[INDEX(cellX, cellY)];

    u[INDEX(cellX, cellY)] = (vx != 0) ? lerp((float) vx, 
      (float) dx, 0.85f) : dx;
    v[INDEX(cellX, cellY)] = (vy != 0) ? lerp((float) vy, 
      (float) dy, 0.85f) : dy;
  }

  void tick(double dt, double visc, double diff) {
    vel_step(u, v, u_prev, v_prev, visc, dt);
  }

  // method used to be 'static' since this class is not a top level type
  final int INDEX(int i, int j) {
    return i + (N + 2) * j;
  }

  double[] tmp = new double[SIZE];

  // same applies to the swap operation ^^
  final void SWAP(double[] x0, double[] x) {
   arraycopy(x0, 0, tmp, 0, SIZE);
   arraycopy(x, 0, x0, 0, SIZE);
   arraycopy(tmp, 0, x, 0, SIZE);
  }

  void add_source(double[] x, double[] s, double dt) {
    int i, size = (N + 2) * (N + 2);
    for (i = 0; i < size; i++)
      x[i] += dt * s[i];
  }

  void diffuse(int b, double[] x, double[] x0, double diff, double dt) {
    int i, j, k;
    double a = dt * diff * N * N;
    for (k = 0; k < 20; k++) {
      for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
          x[INDEX(i, j)] = (x0[INDEX(i, j)] + a
            * (x[INDEX(i - 1, j)] + x[INDEX(i + 1, j)]
            + x[INDEX(i, j - 1)] + x[INDEX(i, j + 1)]))
            / (1 + 4 * a);
        }
      }
      set_bnd(b, x);
    }
  }

  void advect(int b, double[] d, double[] d0, double[] u, double[] v, 
    double dt) {
    int i, j, i0, j0, i1, j1;
    double x, y, s0, t0, s1, t1, dt0;
    dt0 = dt * N;
    for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        x = i - dt0 * u[INDEX(i, j)];
        y = j - dt0 * v[INDEX(i, j)];
        if (x < 0.5)
          x = 0.5;
        if (x > N + 0.5)
          x = N + 0.5;
        i0 = (int) x;
        i1 = i0 + 1;
        if (y < 0.5)
          y = 0.5;
        if (y > N + 0.5)
          y = N + 0.5;
        j0 = (int) y;
        j1 = j0 + 1;
        s1 = x - i0;
        s0 = 1 - s1;
        t1 = y - j0;
        t0 = 1 - t1;
        d[INDEX(i, j)] = s0
          * (t0 * d0[INDEX(i0, j0)] + t1 * d0[INDEX(i0, j1)])
          + s1
          * (t0 * d0[INDEX(i1, j0)] + t1 * d0[INDEX(i1, j1)]);
      }
    }
    set_bnd(b, d);
  }

  void set_bnd(int b, double[] x) {
    int i;
    for (i = 1; i <= N; i++) {
      x[INDEX(0, i)] = (b == 1) ? -x[INDEX(1, i)] : x[INDEX(1, i)];
      x[INDEX(N + 1, i)] = b == 1 ? -x[INDEX(N, i)] : x[INDEX(N, i)];
      x[INDEX(i, 0)] = b == 2 ? -x[INDEX(i, 1)] : x[INDEX(i, 1)];
      x[INDEX(i, N + 1)] = b == 2 ? -x[INDEX(i, N)] : x[INDEX(i, N)];
    }
    x[INDEX(0, 0)] = 0.5 * (x[INDEX(1, 0)] + x[INDEX(0, 1)]);
    x[INDEX(0, N + 1)] = 0.5 * (x[INDEX(1, N + 1)] + x[INDEX(0, N)]);
    x[INDEX(N + 1, 0)] = 0.5 * (x[INDEX(N, 0)] + x[INDEX(N + 1, 1)]);
    x[INDEX(N + 1, N + 1)] = 0.5 * (x[INDEX(N, N + 1)] + x[INDEX(N + 1, N)]);
  }

  void dens_step(double[] x, double[] x0, double[] u, double[] v, 
    double diff, double dt) {
    add_source(x, x0, dt);
    SWAP(x0, x);
    diffuse(0, x, x0, diff, dt);
    SWAP(x0, x);
    advect(0, x, x0, u, v, dt);
  }

  void vel_step(double[] u, double[] v, double[] u0, double[] v0, 
    double visc, double dt) {
    diffuse(1, u, u, visc, dt);
    diffuse(2, v, v, visc, dt);
    project(u, v, u0, v0);
  }

  void project(double[] u, double[] v, double[] p, double[] div) {
    int i, j, k;
    double h;
    h = 1.0 / N;
    for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        div[INDEX(i, j)] = -0.5
          * h
          * (u[INDEX(i + 1, j)] - u[INDEX(i - 1, j)]
          + v[INDEX(i, j + 1)] - v[INDEX(i, j - 1)]);
        p[INDEX(i, j)] = 0;
      }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    for (k = 0; k < 20; k++) {
      for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
          p[INDEX(i, j)] = (div[INDEX(i, j)] + p[INDEX(i - 1, j)]
            + p[INDEX(i + 1, j)] + p[INDEX(i, j - 1)] + p[INDEX(
            i, j + 1)]) / 4;
        }
      }
      set_bnd(0, p);
    }
    for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        u[INDEX(i, j)] -= 0.5
          * (p[INDEX(i + 1, j)] - p[INDEX(i - 1, j)]) / h;
        v[INDEX(i, j)] -= 0.5
          * (p[INDEX(i, j + 1)] - p[INDEX(i, j - 1)]) / h;
      }
    }
    set_bnd(1, u);
    set_bnd(2, v);
  }
}

//Reset coffee and particle colors
void mousePressed(){
	frameCount=0;
}
