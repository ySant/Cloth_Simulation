/* 
 *  Cloth Simulation
 *  --
 *  Yudi Santoso, 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include <GL/glut.h>


// Cloth variables (nv x nh vertex array)
static const GLuint nh = 16;  // horizontal dimension
static const GLuint nv = 16;  // vertical dimension
static GLfloat x[nh][nv][3]; 
static GLfloat v[nh][nv][3]; 
static GLfloat F[nh][nv][3]; 

// The mesh spring parameters:
static const  GLfloat k_str = 3.5f;      // structural
static const  GLfloat k_sh = 0.8f;      // shear
static const  GLfloat k_bnd = 0.2f;    // bending
static const  GLfloat l0h = 0.1f;     // rest length - horizontal
static const  GLfloat l0v = 0.1f;     // rest length - vertical
static const  GLfloat l0d = (GLfloat) sqrt(l0h*l0h+l0v*l0v);   // rest length - diagonal
static const  GLfloat m = 1.0f;   // mass for each point (uniform)
static const  GLfloat d = 0.1f;   // damping constant (uniform for simplicity)

// Wind:
static GLfloat wind_x = -0.0002;
static GLfloat wind_z = 0.0025;
// Gravity
static const GLfloat grav = -0.0012f;

static struct timeval prev_time;   // time 

static const int sPF = 10;  // steps per frame
static float dt = 0.0;      // time interval - to be set later



double drand(double dMin, double dMax)  // generate random double
{
  double a = (double) rand()/RAND_MAX;
  return dMin + a*(dMax-dMin);
}

void fixCorners()
{
// Fix the two upper corners:
      x[0][0][0] = -0.4f + 0.02f;
      x[0][0][1] = 0.9f;
      x[0][0][2] = 0.0f;
      x[nh-1][0][0] = -0.3f + (nh-1)*(1.0f/nh) - 0.02f;
      x[nh-1][0][1] = 0.9f;
      x[nh-1][0][2] = 0.0f;
      v[0][0][0] = 0.0f;
      v[0][0][1] = 0.0f;
      v[0][0][2] = 0.0f;
      v[nh-1][0][0] = 0.0f;
      v[nh-1][0][1] = 0.0f;
      v[nh-1][0][2] = 0.0f;
}

void initialState()  // set initial positions and velocities 
{
// Initial position: - randomize within a range
  double rd;
// random seed:
  srand (time(NULL));
  
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      rd = drand(-0.03, 0.03);
      x[i][j][0] = -0.4f + i*(1.0f/nh); // + rd;
      rd = drand(-0.03, 0.03);
      x[i][j][1] = 0.9f - j*(1.0f/nv); // + rd;
      rd = drand(-0.03, 0.03);
      x[i][j][2] = 0.0f;// + rd;
    }
  }

// Initial velocity: - randomize within a range

  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      v[i][j][0] = drand(-0.001, 0.001);
      v[i][j][1] = drand(-0.001, 0.001);
      v[i][j][2] = drand(-0.001, 0.001);
    }
  }

  fixCorners();

}

void updateForce()
{
// The forces - including the edges effect:

  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        F[i][j][k] = 0;
      }
// add gravity
        F[i][j][1] += grav;
// add wind
    //  if (wind_x < 0) wind_x += 0.00001;
        F[i][j][0] += wind_x;
   //   if (wind_z > 0) wind_z += -0.00001;
        F[i][j][3] += wind_z;
        

      if (i>0) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i-1][j][k]), 2.0);
          dvx += (v[i][j][k]-v[i-1][j][k])*(x[i][j][k]-x[i-1][j][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_str*(x[i][j][k]-x[i-1][j][k])/lt*(lt-l0h);
          F[i][j][k] += -d*(x[i][j][k]-x[i-1][j][k])/dlt*dvx;
        }
      }

      if (j>0) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i][j-1][k]), 2.0);
          dvx += (v[i][j][k]-v[i][j-1][k])*(x[i][j][k]-x[i][j-1][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_str*(x[i][j][k]-x[i][j-1][k])/lt*(lt-l0v);
          F[i][j][k] += -d*(x[i][j][k]-x[i][j-1][k])/dlt*dvx;
        }
      }

      if ((i>0) && (j>0)) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i-1][j-1][k]), 2.0);
          dvx += (v[i][j][k]-v[i-1][j-1][k])*(x[i][j][k]-x[i-1][j-1][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_sh*(x[i][j][k]-x[i-1][j-1][k])/lt*(lt-l0d);
          F[i][j][k] += -d*(x[i][j][k]-x[i-1][j-1][k])/dlt*dvx;
        }
      }

      if (i<nh-1) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i+1][j][k]), 2.0);
          dvx += (v[i][j][k]-v[i+1][j][k])*(x[i][j][k]-x[i+1][j][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_str*(x[i][j][k]-x[i+1][j][k])/lt*(lt-l0h);
          F[i][j][k] += -d*(x[i][j][k]-x[i+1][j][k])/dlt*dvx;
        }
      }

      if (j<nv-1) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i][j+1][k]), 2.0);
          dvx += (v[i][j][k]-v[i][j+1][k])*(x[i][j][k]-x[i][j+1][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_str*(x[i][j][k]-x[i][j+1][k])/lt*(lt-l0v);
          F[i][j][k] += -d*(x[i][j][k]-x[i][j+1][k])/dlt*dvx;
        }
      }

      if ((i<nh-1) && (j<nv-1)) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i+1][j+1][k]), 2.0);
          dvx += (v[i][j][k]-v[i+1][j+1][k])*(x[i][j][k]-x[i+1][j+1][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_sh*(x[i][j][k]-x[i+1][j+1][k])/lt*(lt-l0d);
          F[i][j][k] += -d*(x[i][j][k]-x[i+1][j+1][k])/dlt*dvx;
       }
      }

      if ((i>0) && (j<nv-1)) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i-1][j+1][k]), 2.0);
          dvx += (v[i][j][k]-v[i-1][j+1][k])*(x[i][j][k]-x[i-1][j+1][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_sh*(x[i][j][k]-x[i-1][j+1][k])/lt*(lt-l0d);
          F[i][j][k] += -d*(x[i][j][k]-x[i-1][j+1][k])/dlt*dvx;
        }
      }

      if ((i<nh-1) && (j>0)) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i+1][j-1][k]), 2.0);
          dvx += (v[i][j][k]-v[i+1][j-1][k])*(x[i][j][k]-x[i+1][j-1][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_sh*(x[i][j][k]-x[i+1][j-1][k])/lt*(lt-l0d);
          F[i][j][k] += -d*(x[i][j][k]-x[i+1][j-1][k])/dlt*dvx;
        }
      }

      if (i>1) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i-2][j][k]), 2.0);
          dvx += (v[i][j][k]-v[i-2][j][k])*(x[i][j][k]-x[i-2][j][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_bnd*(x[i][j][k]-x[i-2][j][k])/lt*(lt-2.0*l0h);
          F[i][j][k] += -d*(x[i][j][k]-x[i-2][j][k])/dlt*dvx;
       }
      }

      if (j>1) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i][j-2][k]), 2.0);
          dvx += (v[i][j][k]-v[i][j-2][k])*(x[i][j][k]-x[i][j-2][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_bnd*(x[i][j][k]-x[i][j-2][k])/lt*(lt-2.0*l0v);
          F[i][j][k] += -d*(x[i][j][k]-x[i][j-2][k])/dlt*dvx;
        }
      }

      if (i<nh-2) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i+2][j][k]), 2.0);
          dvx += (v[i][j][k]-v[i+2][j][k])*(x[i][j][k]-x[i+2][j][k]);
        }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_bnd*(x[i][j][k]-x[i+2][j][k])/lt*(lt-2.0*l0h);
          F[i][j][k] += -d*(x[i][j][k]-x[i+2][j][k])/dlt*dvx;
        }
      }

      if (j<nv-2) {
        double dlt = 0;
        double dvx = 0;
        for (int k=0; k<3; k++){
          dlt += pow((x[i][j][k]-x[i][j+2][k]), 2.0);
           dvx += (v[i][j][k]-v[i][j+2][k])*(x[i][j][k]-x[i][j+2][k]);
       }
        GLfloat lt = (GLfloat)(sqrt(dlt));
        for (int k=0; k<3; k++){
          F[i][j][k] += -k_bnd*(x[i][j][k]-x[i][j+2][k])/lt*(lt-2.0*l0v);
          F[i][j][k] += -d*(x[i][j][k]-x[i][j+2][k])/dlt*dvx;
        }
      }

    }
  }

}

void updateState_explicitEuler()
{
// Explicit Euler
  updateForce();

  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        x[i][j][k] = x[i][j][k] + v[i][j][k]*dt;
        v[i][j][k] = v[i][j][k] + F[i][j][k]/m*dt;
      }
    }
  }

}

// Runge-Kutta
void updateState_RK2()
{
// K1:
  GLfloat K1x[nh][nv][3];   // xt
  GLfloat K1v[nh][nv][3];   // vt
  GLfloat K1a[nh][nv][3];   // at
  updateForce();
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K1x[i][j][k] = x[i][j][k];   
        K1v[i][j][k] = v[i][j][k];
        K1a[i][j][k] = F[i][j][k]/m;
      }
    }
  }
// K2:
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        x[i][j][k] = x[i][j][k] + 0.5*dt*K1v[i][j][k];
        v[i][j][k] = v[i][j][k] + 0.5*dt*K1a[i][j][k];  // K2v
      }
    }
  }
  updateForce();   // F at t+dt/2,  K2a = F(t+dt/2)/m
     
// Update state:
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        x[i][j][k] = K1x[i][j][k] + dt*v[i][j][k];
        v[i][j][k] = K1v[i][j][k] + dt*F[i][j][k]/m;
      }
    }
  }

}

void updateState_RK4()
{
// K1:
  GLfloat K1x[nh][nv][3];   // xt
  GLfloat K1v[nh][nv][3];   // vt
  GLfloat K1a[nh][nv][3];   // at
  updateForce();            // Ft
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K1x[i][j][k] = x[i][j][k];   
        K1v[i][j][k] = v[i][j][k];
        K1a[i][j][k] = F[i][j][k]/m;
      }
    }
  }
// K2:
  GLfloat K2x[nh][nv][3];   
  GLfloat K2v[nh][nv][3];   
  GLfloat K2a[nh][nv][3];   
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K2x[i][j][k] = x[i][j][k] + 0.5*dt*K1v[i][j][k];
        K2v[i][j][k] = v[i][j][k] + 0.5*dt*K1a[i][j][k];  
        x[i][j][k] = K2x[i][j][k];   // to compute F
        v[i][j][k] = K2v[i][j][k];   // to compute F
      }
    }
  }
  updateForce();   // F at t+dt/2,  K2a = F(S2, t+dt/2)/m
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K2a[i][j][k] = F[i][j][k]/m;
      }
    }
  }
// K3:
  GLfloat K3x[nh][nv][3];   
  GLfloat K3v[nh][nv][3];   
  GLfloat K3a[nh][nv][3];   
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K3x[i][j][k] = K1x[i][j][k] + 0.5*dt*K2v[i][j][k];
        K3v[i][j][k] = K1v[i][j][k] + 0.5*dt*K2a[i][j][k];  
        x[i][j][k] = K3x[i][j][k];   // to compute F
        v[i][j][k] = K3v[i][j][k];   // to compute F
      }
    }
  }
  updateForce();   // F at t+dt/2,  K3a = F(S3, t+dt/2)/m
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K3a[i][j][k] = F[i][j][k]/m;
      }
    }
  }
// K4:
  GLfloat K4x[nh][nv][3];   
  GLfloat K4v[nh][nv][3];   
  GLfloat K4a[nh][nv][3];   
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K4x[i][j][k] = K1x[i][j][k] + dt*K3v[i][j][k];
        K4v[i][j][k] = K1v[i][j][k] + dt*K3a[i][j][k];  
        x[i][j][k] = K4x[i][j][k];   // to compute F
        v[i][j][k] = K4v[i][j][k];   // to compute F
      }
    }
  }
  updateForce();   // F at t+dt,  K3a = F(S4, t+dt)/m
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        K4a[i][j][k] = F[i][j][k]/m;
      }
    }
  }
// Update state:
  for (int i=0; i<nh; i++){
    for (int j=0; j<nv; j++){
      for (int k=0; k<3; k++){
        x[i][j][k] = K1x[i][j][k] + dt/6*(K1v[i][j][k]
                     +2.0*K2v[i][j][k] + 2.0*K3v[i][j][k]
                     +K4v[i][j][k]);
        v[i][j][k] = K1v[i][j][k] + dt/6*(K1a[i][j][k]
                     +2.0*K2a[i][j][k] + 2.0*K3a[i][j][k]
                     +K4a[i][j][k]);
      }
    }
  }

}


void initGL() {
  glClearColor(0.0f, 0.0f, 0.1f, 1.0f);
}

void drawMesh()
{
  for (int i=0; i<nh-1; i++) {
    for (int j=0; j<nv-1; j++) {
      glBegin(GL_LINE_LOOP);
      glColor3f(0.0f, 1.0f, 0.0f);
      glVertex2f(x[i][j][0], x[i][j][1]);
      glVertex2f(x[i+1][j][0], x[i+1][j][1]);
      glVertex2f(x[i][j+1][0], x[i][j+1][1]);
      glEnd();
      glBegin(GL_LINE_LOOP);
      glColor3f(0.0f, 1.0f, 0.0f);
      glVertex2f(x[i+1][j+1][0], x[i+1][j+1][1]);
      glVertex2f(x[i+1][j][0], x[i+1][j][1]);
      glVertex2f(x[i][j+1][0], x[i][j+1][1]);
      glEnd();
    }
  }
}

void drawPoints()
{
  glBegin(GL_POINTS);
  glColor3f(1.0f, 0.0f, 0.0f);
  
  for (int i=0; i<nh; i++) {
    for (int j=0; j<nv; j++) {
      glVertex2f(x[i][j][0], x[i][j][1]);
    }
  }
  
  glEnd();
}


void display()
{
  glClear(GL_COLOR_BUFFER_BIT);
//  drawPoints();
  drawMesh();
  glFlush();
}

void AnimateScene()
{
// Time elapsed 
  struct timeval current_time;
  gettimeofday(&current_time, NULL);
  dt = (float)(current_time.tv_sec  - prev_time.tv_sec) +
  1.0e-6*(current_time.tv_usec - prev_time.tv_usec);

// Animate the motion of the cloth
  for (int s=0; s<sPF; s++)
  {
// Choose the integrator:
//    updateState_explicitEuler();
//    updateState_RK2();
    updateState_RK4();

// Fix the two upper corner:
    fixCorners();
  }

// Save time_now for next time
  prev_time = current_time;
// Force redraw
  glutPostRedisplay();
}

int main()
{
  int argc = 1;
  int t=0;
  char *argv[1] = {(char*)"None"};
  glutInit(&argc, argv);  
  glutInitWindowSize(600, 600);
  glutInitWindowPosition(50, 50);
  glutCreateWindow("Cloth Simulation");


  glutDisplayFunc(display);
  glutIdleFunc (AnimateScene);
  initGL();
  initialState();
  gettimeofday (&prev_time, NULL);
  glutMainLoop();

  return 0;
}

