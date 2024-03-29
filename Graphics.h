// We use the cairo graphics library
// documentation: https://www.cairographics.org/manual/cairo-Paths.html
#include <unistd.h>
#include <assert.h>
#include <cairo/cairo.h>
#include <cairo/cairo-xlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

typedef struct {
  double x,y,vx,vy;
} Particle; //this is our definition of  a particle, position plus speed


class Graphics{
 private:
  int Np;
  double lmin, lmax,Dim,diam;
  cairo_surface_t *sfc;
  cairo_t *cr;
  Display *dsp;
  Drawable da;
 public:
  Graphics(int N, int Pix, double dmin, double dmax,double diameter){
    Np=N;
    assert( dmin < dmax );
    assert(Np >0) ;
    lmin=dmin;
    lmax=dmax;
    Dim=Pix;
    diam=diameter;

    if (( dsp = XOpenDisplay(NULL)) == NULL)      exit(1);//window management X11, and cairo graphics
    int screen = DefaultScreen(dsp);
    da = XCreateSimpleWindow(dsp, DefaultRootWindow(dsp), 0, 0, Dim, Dim, 0, 0, 0);
    XMapWindow(dsp, da);
    sfc = cairo_xlib_surface_create(dsp, da, DefaultVisual(dsp, screen), Dim , Dim);
    cairo_xlib_surface_set_size(sfc, Dim, Dim);
    cr = cairo_create (sfc);
  } //empty window now on screen

  //some error messages on trying to duplicate the window
  Graphics & operator=(Graphics &g){// stop the program when copying the window
    printf("Don't use = with graphics objects  %p \n", &g ); exit(1);
  }
  Graphics(const Graphics &g ){// stop the program when passing windows as argument
    printf("Don't pass graphics objects: use a pointer  %p \n" , &g); exit(2);
  }
  
  
  void draw(Particle*);//draw the particles
  void frame(double , double , double , double );//draw a square
  ~Graphics(){cairo_destroy (cr);cairo_surface_destroy (sfc); } //clean up function
};
