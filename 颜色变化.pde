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

public void setup() {
  fullScreen();  // 使用全屏模式
  fluidSolver = new NavierStokesSolver();
  numParticles = (int)pow(2, 15);
  particles = new Particle[numParticles];
  visc = 0.00025f;
  diff = 0.03f;
  velocityScale = 10;
  vectors = true;
  limitVelocity = 5;
  initParticles();
  rectMode(CENTER);
  
  background(0);  // 初始黑色背景
  rectMode(CORNER);  // 确保rect从左上角开始绘制
}

private void initParticles() {
  for (int i = 0; i < numParticles - 1; i++) {
    particles[i] = new Particle();
    PVector initialPos = PVector.random2D();
    initialPos.mult(random(height/6));
    particles[i].x = width/2 + int(initialPos.x);
    particles[i].y = height/2 + int(initialPos.y);
  }
}

public void draw() {
	if (frameCount%180 == 0){
	  changeX = random(width);
	  changeY = random(height);
  }
	
	if (mouseX < changeX) mouseX++;
	if (mouseX > changeX) mouseX--;
	if (mouseY < changeY) mouseY++;
	if (mouseY > changeY) mouseY--;
  handleMouseMotion();
  
  // 使用半透明的黑色背景来创建拖尾效果
  fill(0, 20);  // 黑色背景，透明度为20
  rect(0, 0, width, height);  // 覆盖整个画布
  
  // Blend colors
  float dissolve = map(frameCount, 0,2000,0,1);
  dissolve = constrain(dissolve, 0,1);
  c1 = color(255, 0, 255);     // 亮粉色
  c2 = color(65, 105, 225);    // 皇家蓝
  c3 = color(0);               // 黑色背景
  c = lerpColor(c1, c2, dissolve*0.9);  // 粒子颜色在粉色和蓝色之间渐变
  surface = c3;
  fill(surface);
  
  double dt = 1 / frameRate;
  fluidSolver.tick(dt, visc, diff);
  vScale = velocityScale * 60. / frameRate;
  drawParticlesPixels();
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
      float lX = p.x - cellX * cellWidth - cellWidth / 2;
      float lY = p.y - cellY * cellHeight - cellHeight / 2;
      int vCell, hCell, vf, hf;
      if (lX > 0) {
        vCell = Math.min(n, cellX + 1);
        vf = 1;
      } 
      else {
        vCell = Math.max(0, cellX - 1);
        vf = -1;
      }
      if (lY > 0) {
        hCell = Math.min(n, cellY + 1);
        hf = 1;
      } 
      else {
        hCell = Math.max(0, cellY - 1);
        hf = -1;
      }
      float dxv = (float) fluidSolver.getDx(vCell, cellY);
      float dxh = (float) fluidSolver.getDx(cellX, hCell);
      float dxvh = (float) fluidSolver.getDx(vCell, hCell);
      float dyv = (float) fluidSolver.getDy(vCell, cellY);
      float dyh = (float) fluidSolver.getDy(cellX, hCell);
      float dyvh = (float) fluidSolver.getDy(vCell, hCell);
      dx = lerp(lerp(dx, dxv, vf * lX / cellWidth), 
      lerp(dxh, dxvh, vf * lX / cellWidth), 
      hf * lY / cellHeight);
      dy = lerp(lerp(dy, dyv, vf * lX / cellWidth), 
      lerp(dyh, dyvh, vf * lX / cellWidth), 
      hf * lY / cellHeight);
      p.x += dx * vScale;
      p.y += dy * vScale;
			
      // 修改边界检查，让粒子在触及画布边缘时重生
      if (p.x < 0 || p.x > width || p.y < 0 || p.y > height) {
        PVector v = PVector.random2D();
        v.mult(random(height/6));
        p.x = width/2 + int(v.x);
        p.y = height/2 + int(v.y);
      }
      
      // 基于粒子位置计算颜色
      float colorMix = (sin(p.x * 0.01) + cos(p.y * 0.01)) * 0.5 + 0.5;
      color particleColor = lerpColor(c1, c2, colorMix);
      set((int) p.x, (int) p.y, particleColor);
    }
  }
}

private void handleMouseMotion() {
  mouseX = max(1, mouseX);
  mouseY = max(1, mouseY);
  int n = NavierStokesSolver.N;
  float cellHeight = height / n;
  float cellWidth = width / n;
	float mouseDx = (mouseX - pmouseX)/3;
  float mouseDy = (mouseY - pmouseY)/3;
  int cellX = floor(mouseX / cellWidth);
  int cellY = floor(mouseY / cellHeight);
  fluidSolver.applyForce(cellX, cellY, mouseDx, mouseDy);
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
