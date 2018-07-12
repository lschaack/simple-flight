/*
 *  Luke Schaack
 *  Homework 3: Lighting and Textures
 *
 *  Sources: Probably every single example, plus google
 *      old fence object: https://www.turbosquid.com/FullPreview/Index.cfm/ID/733492
 *
 *  Key bindings:
 *  mouse      Look and fly around
 *  arrows     Look and fly around, much slower than the mouse
 *  r/R        Back to starting position
 * 
 *  space      Give the airplane a "push," can be held down to act as a sort of
 *                throttle, but the paper airplane texture will glitch back and
 *                forth while it is held.
 *
 *  qe         Scale physics by orders of 2 (q divide, e multiply)
 *
 *  f/F        Toggle fog
 *  s/S        Toggle between short (18s) and long (360s) day
 *  'l'        Toggle lighting
 *  'L'        Toggle smooth/flat lighting
 * 
 *  1          Look at the traffic light
 *  2          Look at the grass
 *  3          Look at the window (third floor of apartments)
 *  0          Look back at the paper airplane
 * 
 *  ESC        Exit
 */
#include "CSCIx229.h"

#define PLASTER 0
#define WOOD 1
#define ALUMINUM 2
#define ASPHALT 3
#define PLASTIC 4
#define GLASS 5
#define GRASS 6
#define GLOW 7
#define DAY 8
#define NIGHT 9
#define PAPER 10

// For normal vectors
#define IN -1
#define OUT 1

#define min(x, y) x < y ? x : y

int zh=0;           //  Azimuth of light
int fov=55;         //  Field of view (for perspective)
double asp=1;       //  Aspect ratio
int dim=36.0;       //  Size of world
int print = 0;      //  Whether to print important vars in the bottom left

/********** Textures, Objects, and Shaders **********/
unsigned int plaster,siding,shingle,hardwood,tabletop,soda,lid,road,xing,
             pavement,fence,plastic,sod,metallic,paint,slats,door,sky;  //  Textures
unsigned int grass, fence_obj;  //  Objects
int grass_shader;
int n_grass = 1024;             // the max number of grasses in the scene
int grass_rot[1024];            // for storing grass rotation, should be length n_grass
GLint tex_loc;

/********** Lighting **********/
int short_day = 0;
int smooth = 1;
int light = 1;
int fog = 0;
float yFloor = -0.001; // for shadows
double light_pos[3] = {0, 2.5, 0};

/********** Useful colors **********/
const float white[] = {1.00,1.00,1.00};
const float brown[] = {0.62,0.32,0.18};
const float beige[] = {0.95,0.95,0.85};
const float metal[] = {0.95,0.98,1.00};
const float gray[] = {0.6, 0.6, 0.6};
const float dim_metal[] = {0.40,0.45,0.50};
const float concrete[] = {0.72,0.68,0.64};
const float sky_blue[] = {0.53,0.80,0.92};
const float sky_blue_dim[] = {0.53,0.80,0.92, 0.05};

const float std_emis[4] = {0, 0, 0, 1.0};

/********** Flight **********/
int push = 1; // start off with a push, just like you actually fly a paper airplane
// determines local airplane axes
double airplane[4][4] = {{1, 0, 0, 0},      // forward
                         {0, 1, 0, 0},      // up
                         {0, 0, 1, 0},      // right
                         {0, 0, 0, 1}};     // not used, for compatability w/MultMatrix
const int terminal_velocity = 4;
double pos[3] = {0, 15, -15};               // position
int pitch = 0, roll = 0, yaw = 0;         // for calculating airplane[4][4] values
int eff_pitch, eff_roll;                    // for improved physics
int X=0, Y=0;
double j_vec[3] = {0, 1, 0}; // for calculating effective pitch and roll

/********** Physics **********/
const double gravity = 9.8;
// weight is mass of a sheet of letter paper in kg * acceleration from gravity
const double weight = 0.0226796 * 9.8;
// Coefficients from boeing 787
const double lift_coeff = 0.01;
const double drag_coeff = 0.012;
double lift = 0, thrust = 0, drag = 0, glide = 0;
double push_val = 400;
double velocity[3] = {0, 0, 0};
double acceleration[3] = {0, 0, 0};

double time_scale = 1.0; // "slows" the physics engine

double t = 0;
double t_prev;
double t_diff = 0;

/********** Camera **********/
int look_light = 0;
int look_grass = 0;
int look_window = 0;
int ph = 0, th = 90;

void printVector3d(double to_print[3]) {
    printf("{%.2f  %.2f  %.2f}\n", to_print[0], to_print[1], to_print[2]);
}

/*
 * Restricts the value of x by 'by' such that x can never
 *  be outside the interval [-by, by]
 */
double limit(double x, double by) {
    if (fabs(x) > by) {
        int sign = x >= 0? 1 : -1;
        return sign * by;
    } else {
        return x;
    }
}

/*
 * Calculates the Euclidean length of a given vector
 */
double magnitude3d(double *to_mag) {
    return sqrt(pow(to_mag[0], 2) + pow(to_mag[1], 2) + pow(to_mag[2], 2));
}

/*
 * Normalizes a given array in-place
 */
void normalize3d(double to_normalize[3]) {
    int len = 3; // length of input and output arrays
    double norm = magnitude3d(to_normalize);

    int i;
    for (i = 0; i < len; i++) {
        to_normalize[i] = to_normalize[i] / norm;
    }
}

/*
 * Calculates the dot product a . b
 */
double dotProduct3d(double a[3], double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* 
 * Calculates the cross product a x b,
 * Inserts the result into res in-place
 */
void crossProduct3d(double a[3], double b[3], double *res) {
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}

/* 
 * Calculates the projection of vector to_project onto the plane defined by
 * normal vector normal.
 * Inserts the result into res in-place
 * Math per: https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
 *  because I slacked in calc 3...
 */
void projectPlane3d(double to_project[3], double normal[3], double res[3]) {
    double dproduct = dotProduct3d(to_project, normal);
    double mag_sqrd = pow(magnitude3d(normal), 2);
    double coeff = dproduct / mag_sqrd;

    res[0] = to_project[0] - coeff * normal[0];
    res[1] = to_project[1] - coeff * normal[1];
    res[2] = to_project[2] - coeff * normal[2];
}

/* 
 * Calculates the angle in degrees between a and b
 * That's all
 */
int findAngle3d(double a[3], double b[3]) {
    double numer = dotProduct3d(a, b);
    double denom = magnitude3d(a) * magnitude3d(b);
    return (int) (acos(numer / denom) * 180 / 3.1415926); // convert rads to degs
}

/* Word soup
 * Extracts parallel component of of in direction in, placing result in direction
 */
void parallelComponent(double of[3], double in[3], double direction[3]) {
    double rot_coeff = dotProduct3d(of, in) / dotProduct3d(in, in);
    direction[0] = rot_coeff * in[0];
    direction[1] = rot_coeff * in[1];
    direction[2] = rot_coeff * in[2];
}

/*
 * Rotates 'rot' about 'about' modifying 'rot' in-place
 */
void rotateVector3d(int angle, double rot[3], double about[3]) {
    // math per https://math.stackexchange.com/a/1432182/394285
    /* Since the projection of a vector on to itself leaves its magnitude
     * unchanged, the dot product of any vector with itself is the square
     * of that vector's magnitude, so denominator should never be 0. */
    double rot_coeff = dotProduct3d(rot, about) / dotProduct3d(about, about);
    double rot_ll[] = {rot_coeff * about[0], // parallel component
                       rot_coeff * about[1],
                       rot_coeff * about[2]};

    double rot_tt[] = {rot[0] - rot_ll[0], // orthogonal component
                       rot[1] - rot_ll[1],
                       rot[2] - rot_ll[2]};


    double w[] = {0, 0, 0};
    crossProduct3d(about, rot_tt, w);

    double rot_tt_mag = magnitude3d(rot_tt);
    double x1 = Cos(angle) / rot_tt_mag;
    double x2 = magnitude3d(w) == 0 ? 0 : Sin(angle) / magnitude3d(w);

    double rot_tt_angle[] = {rot_tt_mag * (x1 * rot_tt[0] + x2 * w[0]),
                             rot_tt_mag * (x1 * rot_tt[1] + x2 * w[1]),
                             rot_tt_mag * (x1 * rot_tt[2] + x2 * w[2])};

    rot[0] = rot_tt_angle[0] + rot_ll[0];
    rot[1] = rot_tt_angle[1] + rot_ll[1];
    rot[2] = rot_tt_angle[2] + rot_ll[2];
}

/*
 * Converts two-dimensional, row-major rowMajor into a one-dimensional,
 * col-major array, placing that array over colMajor
 */
void flattenToColMajor4d(double rowMajor[4][4], double colMajor[16]) {
    int i, j;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            // i determines row, j is the column offset
            colMajor[4 * i + j] = rowMajor[i][j];
        }
    }
}

/*
 *  Draw vertex in cylindrical coordinates, y and z are reversed per convention
 *  of the y axis as "up" in the starting view settings
 */
static void CylindricalVertex3d(double radius, double azimuth, double elevation) {
   glVertex3d(radius * Cos(azimuth) , elevation , radius * Sin(azimuth));
}

/*
 * A convenience function for setting common material settings in one line
 * based on macro name.
 */
static void set_material(int material) {
    const float zero_shine = 0.0;
    const float some_shine = 0.5;
    const float full_shine = 1.0;

    float emis[4] = {std_emis[0], std_emis[1], std_emis[2], std_emis[3]};

    if (material == PAPER) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, beige);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, beige);
        glMaterialf(GL_FRONT, GL_SHININESS, zero_shine);
    } else if (material == PLASTER) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, beige);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, beige);
        glMaterialf(GL_FRONT, GL_SHININESS, zero_shine);
    } else if (material == PLASTIC) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, gray);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, gray);
        glMaterialf(GL_FRONT, GL_SHININESS, some_shine);
    } else if (material == GLASS) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, gray);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, gray);
        glMaterialf(GL_FRONT, GL_SHININESS, full_shine);
    } else if (material == GLOW) {
        emis[2] = 0;    // set emission color to yellow
        emis[3] = 0.5;  // set actual emission value to 0.1
        glMaterialfv(GL_FRONT, GL_SPECULAR, gray);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, gray);
        glMaterialf(GL_FRONT, GL_SHININESS, full_shine);
    } else if (material == GRASS) {
        float grass_green[3] = {0.5, 1, 0};
        glMaterialfv(GL_FRONT, GL_SPECULAR, grass_green);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, grass_green);
        glMaterialf(GL_FRONT, GL_SHININESS, some_shine);
    } else if (material == WOOD) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, brown);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, brown);
        glMaterialf(GL_FRONT, GL_SHININESS, some_shine);
    } else if (material == ALUMINUM) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, metal);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, dim_metal);
        glMaterialf(GL_FRONT, GL_SHININESS, full_shine);
    } else if (material == ASPHALT) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, concrete);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, concrete);
        glMaterialf(GL_FRONT, GL_SHININESS, some_shine);
    } else if (material == slats) {
        float slats_red[3] = {0.8, 0.25, 0.33};
        glColor3f(1, 1, 1);
        glMaterialfv(GL_FRONT, GL_SPECULAR, slats_red);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, slats_red);
        glMaterialf(GL_FRONT, GL_SHININESS, zero_shine);
    } else if (material == DAY) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, sky_blue);
        // Make light emission proportional to height of the sun
        float emis_val = Sin(zh);
        emis[3] = emis_val > 0? emis_val : 0;
        glMaterialf(GL_FRONT, GL_SHININESS, zero_shine);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, sky_blue);
    } else if (material == NIGHT) {
        glMaterialfv(GL_FRONT, GL_SPECULAR, sky_blue_dim);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, sky_blue_dim);
        glMaterialf(GL_FRONT, GL_SHININESS, zero_shine);
    }
    // Universally set emission to switch off day/night emission for other materials
    glMaterialfv(GL_FRONT,GL_EMISSION,emis);
}

/*
 * Draws a rectangle composed of n_quads_horiz smaller rectangles in the x direction
 * and n_quads_vert smaller rectangles in the y direction, for use in drawing
 * rectangles in an arbitrary direction which will still light well depending
 * on context. 
 */
static void drawRectangle(float width, float height,
                          int n_quads_horiz, int n_quads_vert,
                          float n_tex_horiz, float n_tex_vert,
                          unsigned int texture, unsigned int rotate) {
    float dx = width / n_quads_horiz;
    float dy = height / n_quads_vert;
    int i, j;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);

    if (rotate) {
        // rotate texture so it's in the right direction
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glTranslatef(0.5,0.5,0.0);
        glRotated(90,0,0,1);
        glTranslatef(-0.5,-0.5,0.0);
        glMatrixMode(GL_MODELVIEW);
    }

    glNormal3d(1, 0, 0);

    glBegin(GL_QUADS);
    for (i = 0; i < n_quads_horiz; i++) {
        for (j = 0; j < n_quads_vert; j++) {
            glTexCoord2f(n_tex_horiz *  i      * dx / width,     n_tex_vert *  j      * dy / height);
            glVertex3f(0 , j * dy, i * dx);

            glTexCoord2f(n_tex_horiz *  i      * dx / width,     n_tex_vert * (j + 1) * dy / height);
            glVertex3f(0 , j * dy + dy , i * dx);

            glTexCoord2f(n_tex_horiz * (i + 1) * dx / width,     n_tex_vert * (j + 1) * dy / height);
            glVertex3f(0 , j * dy + dy , i * dx + dx);

            glTexCoord2f(n_tex_horiz * (i + 1) * dx / width,     n_tex_vert *  j      * dy / height);
            glVertex3f(0 , j * dy, i * dx + dx);
        }
    }
    
    glEnd();

    if (rotate) {
        // un-rotate texture so it's in the right direction
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glTranslatef(0.5,0.5,0.0);
        glRotated(0,0,0,1);
        glTranslatef(-0.5,-0.5,0.0);
        glMatrixMode(GL_MODELVIEW);
    }

    glDisable(GL_TEXTURE_2D);
}

/*
 *  Draw a ball
 *     at (x,y,z)
 *     radius r
 */
static void ball(double x,double y,double z,double r)
{
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(r,r,r);
   //  White ball
   glColor3f(1,1,1);
   glutSolidSphere(1.0,16,16);
   //  Undo transofrmations
   glPopMatrix();
}

/*
 * Draws the 3D analog of a rectangle, with dimensions, composite number of quads
 * per face, number of textures per face, particular texture, and boolean opacity
 * all specified by the caller.
 */
static void rectangularPrism(float height, float width, float depth,
                             int n_quads_horiz, int n_quads_vert,
                             float n_tex_horiz, float n_tex_vert,
                             unsigned int texture, unsigned int opaque) {
    if (!opaque) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(0);
    }
    
    /********** Front and back **********/
    glPushMatrix();
    // no rotations or translations required for front
    drawRectangle(width, height,
                  n_quads_horiz, n_quads_vert,
                  n_tex_horiz, n_tex_vert,
                  texture, 0);
    glPopMatrix();

    // back
    glPushMatrix();
    glTranslated(-depth, 0, width);
    glRotatef(180, 0, 1, 0);
    drawRectangle(width, height,
                  n_quads_horiz, n_quads_vert,
                  n_tex_horiz, n_tex_vert,
                  texture, 0);
    glPopMatrix();

    /********** Edges **********/
    // right
    glPushMatrix();
    glTranslated(-depth, 0, 0);
    glRotatef(90, 0, 1, 0);
    drawRectangle(depth, height,
                  n_quads_horiz, n_quads_vert,
                  n_tex_horiz, n_tex_vert,
                  texture, 0);
    glPopMatrix();

    // left 
    glPushMatrix();
    glTranslated(0, 0, width);
    glRotatef(-90, 0, 1, 0);
    drawRectangle(depth, height,
                  n_quads_horiz, n_quads_vert,
                  n_tex_horiz, n_tex_vert,
                  texture, 0);
    glPopMatrix();

    // top 
    glPushMatrix();
    glTranslated(0, height, 0);
    glRotatef(90, 0, 0, 1);
    drawRectangle(width, depth,
                  n_quads_horiz, n_quads_vert,
                  n_tex_horiz, n_tex_vert,
                  texture, 0);
    glPopMatrix();

    // bottom
    glPushMatrix();
    glTranslated(-depth, 0, 0);
    glRotatef(-90, 0, 0, 1);
    drawRectangle(width, depth,
                  n_quads_horiz, n_quads_vert,
                  n_tex_horiz, n_tex_vert,
                  texture, 0);
    glPopMatrix();

    if (!opaque)
   {
      glDisable(GL_BLEND);
      glDepthMask(1);
   }
}

/*
 * Draws a transparent window
 */
static void window() {
    float depth = 0.15;         // approx. 6 inches
    float thick = 0.05;         //         2 inches
    float width = 1;
    float height = 1.5; // -thick to avoid overlap w/wall

    float half_height = 0.5 * height;
    float half_width = 0.5 * width;

    set_material(PLASTER);
    glColor3f(0.8, 0.8, 0.8);

    // top board
    glPushMatrix();
    glTranslated(0, half_height - thick, 0);
    rectangularPrism(thick, width, depth, 3, 2, 1, 0.2, paint, 1);
    glPopMatrix();
    // bottom board
    glPushMatrix();
    glTranslated(depth, -half_height - thick, 0);
    rectangularPrism(thick, width, depth * 2, 3, 2, 1, 0.2, paint, 1);
    glPopMatrix();
    // right board
    glPushMatrix();
    glTranslated(0, -half_height - thick, width - thick);
    rectangularPrism(height, thick, depth, 2, 3, 0.2, 1, paint, 1);
    glPopMatrix();
    // left board
    glPushMatrix();
    glTranslated(0, -half_height - thick, 0);
    rectangularPrism(height, thick, depth, 2, 3, 0.2, 1, paint, 1);
    glPopMatrix();

    // Glass panel
    if (zh < 180)
        set_material(GLASS);
    else
        set_material(GLOW);
    glColor4f(0, 0.8, 0.4, 0.3);

    glPushMatrix();
    glTranslated(-depth / 2, -half_height - thick, 0);
    rectangularPrism(height, width, thick / 4, 3, 6, 1, 1, paint, 0);
    glPopMatrix();

    set_material(PLASTER);
    glColor3f(0.8, 0.8, 0.8);

    // vertical bar
    glPushMatrix();
    glTranslated(-depth / 4, -half_height - thick, half_width);
    rectangularPrism(height, thick / 2, depth / 2, 2, 3, 0.2, 1, paint, 1);
    glPopMatrix();
    // horizontal bar
    glPushMatrix();
    glTranslated(-depth / 4, 0, 0);
    rectangularPrism(thick / 2, width, depth / 2, 3, 2, 1, 0.2, paint, 1);
    glPopMatrix();
}

/*
 * Makes a wall out of several quads to improve lighting resolution,
 */
static void wall(double in_out, int material, unsigned int texture, unsigned int addWindow) {
    // n_quads is number of quads that will make up the wall ***in each direction***
    /* in_out determines the normal vector:
     * -1 = in
     *  1 = out 
     */
    const float width = 10.0;
    const float height = width / 2;
    // the below settings naturally yield a wall half as high as it is wide
    const int n_quads = 10;
    float h_step = 1; // horizontal step
    float v_step = 0.5; // vertical step

    glEnable(GL_TEXTURE_2D);
    set_material(material);
    glBindTexture(GL_TEXTURE_2D, texture);

    int i, j;
    for (i = 0; i < n_quads; i++) {
        for (j = 0; j < n_quads; j++) {
            // only draw polygon if: addWindow, the wall faces out, and the quad will not take up window space
            if (!addWindow || (i != 2 && i != 7) || (j < 4 || j > 6)) {
                glNormal3f(in_out, 0, 0);
                glBegin(GL_QUADS);
                if (in_out == IN) {
                    glTexCoord2f(0   , 0   ); glVertex3f(0 , j * v_step          , -width + i * h_step);
                    glTexCoord2f(0.25, 0   ); glVertex3f(0 , j * v_step          , -width + i * h_step + h_step);
                    glTexCoord2f(0.25, 0.25); glVertex3f(0 , j * v_step + v_step , -width + i * h_step + h_step);
                    glTexCoord2f(0   , 0.25); glVertex3f(0 , j * v_step + v_step , -width + i * h_step);
                } else {
                    glTexCoord2f(0   , 0   ); glVertex3f(0 , j * v_step          , -width + i * h_step);
                    glTexCoord2f(0   , 0.25); glVertex3f(0 , j * v_step + v_step , -width + i * h_step);
                    glTexCoord2f(0.25, 0.25); glVertex3f(0 , j * v_step + v_step , -width + i * h_step + h_step);
                    glTexCoord2f(0.25, 0   ); glVertex3f(0 , j * v_step          , -width + i * h_step + h_step);
                }
                glEnd();
            }
        }
    }

    glDisable(GL_TEXTURE_2D);

    // leave spaces above, but only add windows if out-facing, or else there'll be two overlapping
    if (addWindow && in_out == OUT) {
        glPushMatrix();
        glTranslated(0, height / 2 + 0.25, -3 * width / 10);
        window();
        glTranslated(0, 0, -5 * width / 10);
        window();
        glPopMatrix();
    }
}

/* 
 * Draw 6m-wide walls, textured with plaster.
 */
static void walls(int in_out, int material, unsigned int texture) {
    int i;
    double width = 10.0;

    glColor3f(0.65, 0.65, 0.55);

    for (i = 0; i < 4; i++) {
        if (i != 1) { // leave one face open for the window
            glPushMatrix();
            glRotatef(90 * i, 0, 1, 0);
            glTranslated(width / 2, 0, width / 2);
            wall(in_out, material, texture, 0);
            glPopMatrix();
        }
    }

    glPushMatrix();
    glRotatef(90, 0, 1, 0);
    glTranslated(width / 2, 0, width / 2);
    wall(in_out, material, texture, 1);
    glPopMatrix();
}

/* 
 * Draw a floor which defines the full size (dim) of the world, textured with
 * wood. Name avoids conflict w/ math.floor.
 */
static void floorTexture(double size, int n_quads, int material, unsigned int texture) {
    // n_quads is the # of quads that will make up the floor ***in each direction***
    const int floor_width = 2 * size;
    float step = floor_width / (float) n_quads;

    set_material(material);

    glEnable(GL_TEXTURE_2D);
    // allow textures to repeat
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glBindTexture(GL_TEXTURE_2D, texture);

    int i, j;
    for (i = 0; i < n_quads; i++) {
        for (j = 0; j < n_quads; j++) {
            glNormal3f(0,1,0);
            glBegin(GL_QUADS);
            glTexCoord2f(1, 0); glVertex3f(-size + i * step , 0 , -size + j * step);
            glTexCoord2f(1, 1); glVertex3f(-size + i * step , 0 , -size + j * step + step);
            glTexCoord2f(0, 1); glVertex3f(-size + i * step + step, 0 , -size + j * step + step);
            glTexCoord2f(0, 0); glVertex3f(-size + i * step + step, 0 , -size + j * step);
            glEnd();
        }
    }
    glDisable(GL_TEXTURE_2D);
}

/* 
 * Draw a ceiling which defines the full size (dim) of the world, textured with plaster.
 */
static void ceiling(double size) {
   const int ceil_width = 2 * size;
   int n_quads = 4; // number of quads that will make up the floor ***in each direction***
   float step = ceil_width / (float) n_quads;

   glEnable(GL_TEXTURE_2D);
   set_material(PLASTER);
   // allow textures to repeat
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   glBindTexture(GL_TEXTURE_2D, plaster);

   int i, j;
   for (i = 0; i < n_quads; i++) {
      for (j = 0; j < n_quads; j++) {
         glNormal3f(0,-1,0);
         glBegin(GL_QUADS);
         glTexCoord2f(0, 0); glVertex3f(-size + i * step + step, size , -size + j * step);
         glTexCoord2f(0, 1); glVertex3f(-size + i * step + step, size , -size + j * step + step);
         glTexCoord2f(1, 1); glVertex3f(-size + i * step , size , -size + j * step + step);
         glTexCoord2f(1, 0); glVertex3f(-size + i * step , size , -size + j * step);
         glEnd();
      }
   }

   glDisable(GL_TEXTURE_2D);
}

/*
 *  Draw dinner table
 *    at (x,y,z)
 */
static void tableLeg(double x, double y, double z, double degrees, double hgt, double wid) {
   // Dimensions used to size table
   const double lwt = wid / 20;   // half cross-sectional _leg _width of _top of table leg
   const double lwb = lwt * 0.6;  // half cross-sectional _leg _width of _bottom of table leg

   glEnable(GL_TEXTURE_2D);
   set_material(WOOD);
   glBindTexture(GL_TEXTURE_2D, tabletop);

   glPushMatrix();
   // rotate/translate are reversed b/c actually yields nice behavior for a square table
   glTranslated(x, y, z);
   glRotatef(degrees, 0, 1, 0); // parameterize rotation so this can do every leg

   glColor3f(beige[0],beige[1],beige[2]);

   glBegin(GL_QUADS);
   // bottom
   glNormal3f(0,-1,0);
   glTexCoord2f(1, 0); glVertex3d( 0   , 0   , 0   );
   glTexCoord2f(1, 1); glVertex3d( lwb , 0   , 0   );
   glTexCoord2f(0, 1); glVertex3d( lwb , 0   , lwb );
   glTexCoord2f(0, 0); glVertex3d( 0   , 0   , lwb );
   // sides
   glNormal3f(-1,0,0);
   glTexCoord2f(0, 0); glVertex3d( 0   , 0   , 0   );
   glTexCoord2f(1, 0); glVertex3d( 0   , 0   , lwb );
   glTexCoord2f(1, 1); glVertex3d( 0   , hgt , lwt );
   glTexCoord2f(0, 1); glVertex3d( 0   , hgt , 0   );
   // -----
   glNormal3f(0,0,1);
   glTexCoord2f(0, 0); glVertex3d( 0   , 0   , lwb );
   glTexCoord2f(1, 0); glVertex3d( lwb , 0   , lwb );
   glTexCoord2f(1, 1); glVertex3d( lwt , hgt , lwt );
   glTexCoord2f(0, 1); glVertex3d( 0   , hgt , lwt );
   // -----
   glNormal3f(1,0,0);
   glTexCoord2f(0, 0); glVertex3d( lwb , 0   , lwb );
   glTexCoord2f(1, 0); glVertex3d( lwb , 0   , 0   );
   glTexCoord2f(1, 1); glVertex3d( lwt , hgt , 0   );
   glTexCoord2f(0, 1); glVertex3d( lwt , hgt , lwt );
   // -----
   glNormal3f(0,0,-1);
   glTexCoord2f(0, 0); glVertex3d( lwb , 0   , 0   );
   glTexCoord2f(1, 0); glVertex3d( 0   , 0   , 0   );
   glTexCoord2f(1, 1); glVertex3d( 0   , hgt , 0   );
   glTexCoord2f(0, 1); glVertex3d( lwt , hgt , 0   );
   
   glEnd();

   //  Undo transformations
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}

/*
 *  Draw dinner table
 *    at (x,y,z)
 *    front towards (dx,dy,dz)
 *    up towards (ux,uy,uz)
 */
static void table(double x, double y, double z, double scale)
{
   // Dimensions used to size table
   const double hgt = 1.0; // _height
   const double wid = 1.5; // _width

   const double thk = 0.1; // _thickness of the top surface
   const double ovg = 0.5; // how much the surface _overhangs the legs

   glEnable(GL_TEXTURE_2D);
   set_material(WOOD);
   // allow textures to repeat
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   glBindTexture(GL_TEXTURE_2D, tabletop);

   glPushMatrix();
   glTranslated(x, y, z);
   // you can't rotate because I said so
   glScaled(scale, scale, scale);

   glColor3f(beige[0],beige[1],beige[2]);
   //  surface
   glBegin(GL_QUADS);
   // very top
   glNormal3f(0,1,0);
   glTexCoord2f(2, 0); glVertex3d( wid + ovg , hgt + thk ,-wid - ovg ); 
   glTexCoord2f(0, 0); glVertex3d(-wid - ovg , hgt + thk ,-wid - ovg );
   glTexCoord2f(0, 2); glVertex3d(-wid - ovg , hgt + thk , wid + ovg );
   glTexCoord2f(2, 2); glVertex3d( wid + ovg , hgt + thk , wid + ovg );
   // bottom of top surface
   glNormal3f(0,-1,0);
   glTexCoord2f(2, 0); glVertex3d( wid + ovg , hgt       , wid + ovg );
   glTexCoord2f(0, 0); glVertex3d(-wid - ovg , hgt       , wid + ovg );
   glTexCoord2f(0, 2); glVertex3d(-wid - ovg , hgt       ,-wid - ovg );
   glTexCoord2f(2, 2); glVertex3d( wid + ovg , hgt       ,-wid - ovg );
   // sides for watertightness, for some reason no iterative approach was working...
   glNormal3f(0,0,1);
   glTexCoord2f(1, 0); glVertex3d( wid + ovg , hgt + thk , wid + ovg );
   glTexCoord2f(1, 1); glVertex3d(-wid - ovg , hgt + thk , wid + ovg );
   glTexCoord2f(0, 1); glVertex3d(-wid - ovg , hgt       , wid + ovg );
   glTexCoord2f(0, 0); glVertex3d( wid + ovg , hgt       , wid + ovg );
   // -----
   glNormal3f(1,0,0);
   glTexCoord2f(1, 0); glVertex3d( wid + ovg , hgt       , wid + ovg );
   glTexCoord2f(1, 1); glVertex3d( wid + ovg , hgt       ,-wid - ovg );
   glTexCoord2f(0, 1); glVertex3d( wid + ovg , hgt + thk ,-wid - ovg );
   glTexCoord2f(0, 0); glVertex3d( wid + ovg , hgt + thk , wid + ovg );
   // -----
   glNormal3f(0,0,-1);
   glTexCoord2f(1, 0); glVertex3d( wid + ovg , hgt       ,-wid - ovg );
   glTexCoord2f(1, 1); glVertex3d(-wid - ovg , hgt       ,-wid - ovg );
   glTexCoord2f(0, 1); glVertex3d(-wid - ovg , hgt + thk ,-wid - ovg );
   glTexCoord2f(0, 0); glVertex3d( wid + ovg , hgt + thk ,-wid - ovg );
   // -----
   glNormal3f(-1,0,0);
   glTexCoord2f(1, 0); glVertex3d(-wid - ovg , hgt       ,-wid - ovg );
   glTexCoord2f(1, 1); glVertex3d(-wid - ovg , hgt       , wid + ovg );
   glTexCoord2f(0, 1); glVertex3d(-wid - ovg , hgt + thk , wid + ovg );
   glTexCoord2f(0, 0); glVertex3d(-wid - ovg , hgt + thk ,-wid - ovg );
   glEnd();

   // legs
   tableLeg(wid, 0, wid, 0, hgt, wid);
   tableLeg(wid, 0, -wid, 90, hgt, wid);
   tableLeg(-wid, 0, -wid, 180, hgt, wid);
   tableLeg(-wid, 0, wid, 270, hgt, wid);

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}

/*
 * Draws a soda can.
 */
static void sodaCan() {
    /*
        * The below two constants determine the number of iterations to perform
        * per loop, and so can be used to dynamically adjust the effective
        * resolution of the can.
    */
    const int d_theta = 36; // The change in theta per iteration, 36 yields 10 iters
    const int n_iter_h = 3; // The absolute number of height iterations

    const double height = 0.9;
    const double top_height = 0.1;
    const double r = 0.3;
    const double dh_bottom = 0.02;
    const double lip_width = 0.05; // how far in (towards center) the top & bottom curves go
    int h_incrementer; // for an int loop that will actually terminate
    double h;          // for doing float math w/value of h_incrementer
    double dh;
    int theta;

    set_material(ALUMINUM);   
    glColor3f(metal[0], metal[1], metal[2]);

    //  Main body
    glBindTexture(GL_TEXTURE_2D, soda);
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUAD_STRIP);
    for (theta = 0; theta <= 360; theta += d_theta) {
        glNormal3f(Cos(theta),
                    0,
                    Sin(theta));

        glTexCoord2f(theta / 360.0f,0); CylindricalVertex3d(r, theta, 0.1);
        glTexCoord2f(theta / 360.0f,1); CylindricalVertex3d(r, theta, height);
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);

    // Top curve
    glBegin(GL_QUAD_STRIP);
    for (h_incrementer = 0; h_incrementer < n_iter_h; h_incrementer++) {
        for (theta = 0; theta <= 360; theta += d_theta) {
            dh = top_height / n_iter_h;
            h = h_incrementer * dh; // h describes the height from the top of the can

            glNormal3f(Cos(theta),
                        Cos((h / (float) 5) * 90),
                        Sin(theta));

            // cos(h*10*90) to get almost 's'-shaped curve
            // 0.1 - h in the lines below to flip the curve vertically
            CylindricalVertex3d(r - (lip_width * Cos(((h + dh) * 900))), theta, height + top_height - h - dh);
            CylindricalVertex3d(r - (lip_width * Cos((h * 900))), theta, height + top_height - h);
        }
    }
    glEnd();

    // lid
    glNormal3f(0,1,0);
    glBegin(GL_QUAD_STRIP);
    // top lip thickness
    for (theta = 0; theta <= 360; theta += d_theta) {
        CylindricalVertex3d(r - lip_width, theta, height + 0.1);
        CylindricalVertex3d(r - lip_width - 0.02, theta, height + 0.1);
    }
    glEnd();
    // top lip inner curve (just the edges of a cylinder)
    glBegin(GL_QUAD_STRIP);
    for (theta = 0; theta <= 360; theta += d_theta) {
        glNormal3f(-Cos(theta), 0, -Sin(theta));
        CylindricalVertex3d(r - lip_width - 0.02, theta, height + 0.1);
        CylindricalVertex3d(r - lip_width - 0.02, theta, height + 0.05);
    }
    glEnd();

    glBindTexture(GL_TEXTURE_2D, lid);
    glEnable(GL_TEXTURE_2D);
    glNormal3f(0,1,0);

    glTexCoord2f(0.5,0.5);
    CylindricalVertex3d(0, 0, height + 0.05);

    glBegin(GL_TRIANGLE_FAN);
    for (theta = 360; theta >= 0; theta -= d_theta) {
        glTexCoord2f(0.5 * Sin(theta) + 0.5, 0.5 * Cos(theta) + 0.5);
        CylindricalVertex3d(r - lip_width - 0.02, theta + d_theta, height + 0.05);
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);

    // Bottom curve
    glBegin(GL_QUAD_STRIP);
    for (h_incrementer = 0; h_incrementer <= n_iter_h; h_incrementer++) {
        for (theta = 0; theta <= 360; theta += d_theta) {
            // h describes the height from the bottom of the can, so range should be [0, 0.1]
            dh = top_height / n_iter_h;
            h = h_incrementer * dh; // h describes the height from the top of the can
            // circular normal w/y component defined by how far we are "into" the range of the Cos curve
            glNormal3f(Cos(theta),
                    -Cos((h / (float) n_iter_h) * 90),
                        Sin(theta));

            // cos(h*n_iter_h*90) to get almost 's'-shaped curve
            // 0.1 - h in the lines below to flip the curve vertically
            CylindricalVertex3d(r - lip_width + (lip_width * Cos((h * 900))), theta, 0.1 - h);
            CylindricalVertex3d(r - lip_width + (lip_width * Cos(((h - dh) * 900))), theta, 0.1 - h + dh);
        }
    }
    glEnd();
    // bottom lip thickness
    glNormal3f(0,-1,0);
    glBegin(GL_QUAD_STRIP);
    for (theta = 360; theta >= 0; theta -= d_theta) {
        CylindricalVertex3d(r - lip_width, theta, 0.0);
        CylindricalVertex3d(r - lip_width - dh_bottom, theta, 0.0);
    }
    glEnd();
    // hemisphere that makes up the very bottom
    glBegin(GL_QUAD_STRIP);
    for (h_incrementer = 0; h_incrementer <= n_iter_h; h_incrementer++) {
        for (theta = 360; theta >= 0; theta -= d_theta) {
            dh = (r - lip_width) / n_iter_h;
            h = h_incrementer * dh; // h describes the height from the top of the can

            glNormal3f(-Sin(theta) * Sin((h / n_iter_h) * 180),
                        -Cos(theta) * Sin((h / n_iter_h) * 180),
                                    -Sin((h / n_iter_h) * 180));

            // sphere in cylindrical coords is R^2 = r^2 + z^2
            CylindricalVertex3d(sqrt(pow(r - lip_width, 2) - pow(h, 2)), theta, h);
            CylindricalVertex3d(sqrt(pow(r - lip_width, 2) - pow(h + dh, 2)), theta, h + dh);
        }
    }
    glEnd();

    glNormal3f(0,-1,0);
    glBegin(GL_TRIANGLE_FAN);
    CylindricalVertex3d(0, 0, r - lip_width);
        for (theta = 0; theta <= 360; theta += d_theta) {
            CylindricalVertex3d(sqrt(pow(r - lip_width, 2) - pow(r-lip_width, 2)), theta, r - lip_width);
            CylindricalVertex3d(sqrt(pow(r - lip_width, 2) - pow(r-lip_width, 2)), theta + d_theta, r - lip_width);
        }
    glEnd();
}

/*
 * Draws a room with pre-positioned table and soda cans
 */
static void drawRoom(int drawFurniture) {
    float room_size = 5.0;

    glPushMatrix();
    glTranslated(0, 0, 0);
    walls(IN, PLASTER, plaster);
    walls(OUT, WOOD, siding);
    floorTexture(room_size, 4, WOOD, hardwood);
    ceiling(room_size);
    if (drawFurniture)
        table(2, 0, 2, 1);
    glPopMatrix();

    if (drawFurniture) {
        glPushMatrix();
        glTranslated(1.3, 1.13, 1.5);
        glScaled(0.3,0.3,0.3);
        sodaCan();
        glPopMatrix();

        glPushMatrix();
        glTranslated(1.2, 1.2, 2);
        glRotatef(90, 1, 0, 0);
        glRotatef(36, 0, 0, 1);
        glScaled(0.3,0.3,0.3);
        sodaCan();
        glPopMatrix();
    }
}


/*
 * Draw the entire apartment by iteratively drawing rooms for both a horizontal
 * and a vertical loop, as well as a roof and wood slat seperators
 */
static void drawApartment() {
    int i, j;
    // Draw 5 columns of 3 rooms high, with slats in-between
    for (i = 0; i < 5; i++) {
        glPushMatrix();
        for (j = 0; j < 3; j++) {
            // only draw furniture in half of rooms for a slight speed increase
            drawRoom(i % 2 != j % 2? 1 : 0);
            glTranslated(0, 5.0, 0);
        }
        glPopMatrix();


        // draw the roof
        // below, all the .99s are done so as not to overlap the siding texture
        glPushMatrix();
        glTranslated(0, 20, -4.99);
        glRotated(45, 0, 0, 1);
        rectangularPrism(0.05, 9.98, sqrt(2) * 5.0 + 1, 3, 6, 6, 4, shingle, 1);
        glPopMatrix();
        // draw the other side
        glPushMatrix();
        glTranslated(5.7, 14.3, -4.99);
        glRotated(-45, 0, 0, 1);
        rectangularPrism(0.05, 9.98, sqrt(2) * 5.0 + 1, 3, 6, 6, 4, shingle, 1);
        glPopMatrix();
        

        // draw the roof face
        glColor3f(0.6, 0.6, 0.5);
        set_material(PLASTER);
        glPushMatrix();
        glTranslated(0, 15, -5);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, siding);
        glBegin(GL_TRIANGLE_FAN);

        glNormal3f(0, 0, -1);
        glTexCoord2d(1, 1); glVertex3d(0, 2.5, 0);

        glTexCoord2d(0, 0); glVertex3d( 5.0, 0, 0);
        glTexCoord2d(2, 0); glVertex3d(-5.0, 0, 0);
        glTexCoord2d(2, 2); glVertex3d(0, 5.0, 0.3);
        glTexCoord2d(0, 0); glVertex3d( 5.0, 0, 0);

        glEnd();
        // other side
        glTranslated(0, 0, 10);
        glBegin(GL_TRIANGLE_FAN);

        glNormal3f(0, 0, 1);
        glTexCoord2d(1, 1); glVertex3d(0, 2.5, 0);

        glTexCoord2d(0, 0); glVertex3d( 5.0, 0, 0);
        glTexCoord2d(2, 2); glVertex3d(   0, 5.0, -0.3);
        glTexCoord2d(2, 0); glVertex3d(-5.0, 0, 0);
        glTexCoord2d(0, 0); glVertex3d( 5.0, 0, 0);

        glEnd();
        glDisable(GL_TEXTURE_2D);
        glPopMatrix();


        if (i < 4) { // don't do this the last loop (none on corners)
            set_material(WOOD);
            glPushMatrix();
            glTranslated(5.0, 0.0, -5.0);
            glRotatef(90, 0, 1, 0);
            drawRectangle(1.0, 15.0, // entire height of 3 rooms, plus space between
                          2, 9, 0.33, 8.5,
                          slats, 0);

            // again on the other side
            glRotatef(180, 0, 1, 0);
            glTranslated(10.0, 0, -1.0);
            drawRectangle(1.0, 15.0,
                          2, 9, 0.33, 8.5,
                          slats, 0);
            glPopMatrix();
        }

        glTranslated(11.0, 0, 0);
    }
}

/*
 * Behaves as strictly a helper function for stopLight()––draws a single light
 * with plastic enclosure
 */
static void drawLight(double x, double y, double z, double r, float color[3], unsigned int light, int on) {
    /********** Draw enclosure **********/
    int theta = 0;
    int d_theta = 20;
    double thick = 0.3;
    double thick2 = thick / 6;
    
    glColor3f(0.3, 0.3, 0.3);
    set_material(PLASTIC);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, PLASTIC);

    glDisable(GL_CULL_FACE); // necessary for viewing 360
    glBegin(GL_QUAD_STRIP);
    for (theta = -60; theta <= 240; theta += d_theta) {
        // theta + 15 to get angle at center of quad
        glNormal3f(0, -Sin(theta), -Cos(theta + 10));

        glTexCoord2f(Sin(theta), 1);  glVertex3d(x + thick, y + r * Sin(theta), z + r * (-Cos(theta)));
        glTexCoord2f(Sin(theta), 0);  glVertex3d(x        , y + r * Sin(theta), z + r * (-Cos(theta)));
    }
    glEnd();
    glEnable(GL_CULL_FACE); // necessary for viewing 360

    glBegin(GL_QUAD_STRIP);
    for (theta = 180; theta <= 360; theta += d_theta) {
        glNormal3f(0, -Sin(theta), -Cos(theta + 10));

        glTexCoord2f(Sin(theta), 1);  glVertex3d(x + thick2, y + r * Sin(theta), z + r * (-Cos(theta)));
        glTexCoord2f(Sin(theta), 0);  glVertex3d(x         , y + r * Sin(theta), z + r * (-Cos(theta)));
    }
    glEnd();

    glDisable(GL_TEXTURE_2D);

    /********** Create the lit panel and light itself **********/
    float position[4] = {x, y, z, 1.0};
    float xDirection[3] = {1, 0, 0};
    float ambient[3] = {0.2, 0.2, 0.2};
    float emission[4] = {color[0], color[1], color[2], 0.1};

    set_material(GLASS);
    // make the panel light up
    glPushMatrix();
    glTranslated(x, y, z);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    if (on) {
        glColor3fv(color);

        glEnable(GL_LIGHT0 + light);

        glLightfv(GL_LIGHT0 + light, GL_POSITION, position);
        glLightfv(GL_LIGHT0 + light, GL_SPOT_DIRECTION, xDirection);

        glLightf(GL_LIGHT0 + light, GL_SPOT_CUTOFF, 45.0);
        glLightf(GL_LIGHT0 + light, GL_LINEAR_ATTENUATION, 0.5);
        glLightf(GL_LIGHT0 + light, GL_SPOT_EXPONENT, 1.0);

        glLightfv(GL_LIGHT0 + light, GL_AMBIENT , ambient);
        glLightfv(GL_LIGHT0 + light, GL_DIFFUSE , color);
        glLightfv(GL_LIGHT0 + light, GL_SPECULAR, color);

        //  Location of viewer for specular calculations
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,1);

        glMaterialfv(GL_FRONT,GL_EMISSION,emission);
        glMaterialfv(GL_FRONT,GL_SPECULAR,color);
    } else {
        glColor3fv(gray);
    }

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(1, 0, 0);
    glVertex3d(0.08, 0, 0);
    for (theta = 0; theta <= 360; theta += d_theta) {
        glVertex3d(thick2, r * Sin(theta), r * (-Cos(theta)));
    }
    glEnd();
    glPopMatrix();
}

/*
 * Draw a stop light with a timed toggle between RYG
 * Behaves non-strictly as a helper function for trafficLight()
 */
static void stopLight(double x, double y, double z, double degrees, double scale) {
    const int height = 1;
    const double width = 0.4;
    const double radius = 0.12;
    const double thick = 0.03;

    glColor3f(0.3, 0.3, 0.3);

    glPushMatrix();
    glTranslated(x, y, z);
    glRotatef(degrees, 0, 1, 0);
    glScaled(scale, scale, scale);

    drawRectangle(width, height, 2, 3, 2, 3, plastic, 0);

    // perform one more matrix push/pop to get backwards orientation
    glPushMatrix();
    glTranslated(-thick, 0, width);
    glRotatef(180, 0, 1, 0);
    drawRectangle(width, height, 2, 3, 2, 3, plastic, 0);
    glPopMatrix();

    /********** edges **********/
    // right
    glPushMatrix();
    glTranslated(-thick, 0, 0);
    glRotatef(90, 0, 1, 0);
    drawRectangle(thick, height, 1, 3, 1, 3, plastic, 0);
    glPopMatrix();

    // left 
    glPushMatrix();
    glTranslated(0, 0, width);
    glRotatef(270, 0, 1, 0);
    drawRectangle(thick, height, 1, 3, 1, 3, plastic, 0);
    glPopMatrix();

    // top 
    glPushMatrix();
    glTranslated(0, height, 0);
    glRotatef(90, 0, 0, 1);
    drawRectangle(width, thick, 3, 1, 3, 1, plastic, 0);
    glPopMatrix();

    // bottom
    glPushMatrix();
    glTranslated(-thick, 0, 0);
    glRotatef(-90, 0, 0, 1);
    drawRectangle(width, thick, 3, 1, 3, 1, plastic, 0);
    glPopMatrix();

    float red[3] = {1, 0, 0};
    float yellow[3] = {1, 1, 0};
    float green[3] = {0, 1, 0};

    // do timing for the automatic lights switch
    int draw_red, draw_yellow, draw_green;
    draw_red = draw_yellow = draw_green = 0;

    float t_fifteen = fmod(t, 15.0);
    if (t_fifteen < 7) {
        draw_green = 1;
    } else if (t_fifteen < 10) {
        draw_yellow = 1;
    } else {
        draw_red = 1;
    }

    drawLight(0, 0.8, width / 2, radius, red, 1, draw_red);
    drawLight(0, 0.5, width / 2, radius, yellow, 2, draw_yellow);
    drawLight(0, 0.2, width / 2, radius, green, 3, draw_green);
    glPopMatrix();
}

/*
 * Draw an entire traffic light––a metal pole with a stoplight on the end
 */
static void trafficLight() {
    set_material(ALUMINUM);
    int i;
    int n_lengths = 10;
    float length_step = 0.5;
    int th;
    int d_th = 36;
    float r = 0.2;

    /********** Vertical pole **********/
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, metallic);
    glBegin(GL_QUAD_STRIP);
    for (i = 0; i < n_lengths; i++) {
        for (th = 360; th >= 0; th -= d_th) {
            glNormal3f(Cos(th), 0, Sin(th));
            glTexCoord2d(Sin(th), 1);  glVertex3d(r * Cos(th), i * length_step + length_step, r * Sin(th));
            glTexCoord2d(Sin(th), 0);  glVertex3d(r * Cos(th), i * length_step              , r * Sin(th));
        }
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);

    // end cap
    glPushMatrix();
    glTranslated(0, n_lengths * length_step, 0);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, metallic);
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 1, 0);
    glTexCoord2d(0.5, 0.5);
    glVertex3f(0, 0.1, 0);
    for (th = 360; th >= 0; th -= d_th) {
        // add a bit of outwardness to the edge lighting so it looks a little less flat
        glNormal3f(Cos(th), 1, Sin(th));
        glTexCoord2d(Cos(th), Sin(th));  glVertex3f(r * Cos(th), 0, r * Sin(th));
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);

    // next section works nicely with this matrix, so no call to glPopMatrix()

    /********** Horizontal crossbar **********/
    r /= 2; // half radius for the crossbar
    n_lengths = 14;
    glTranslated(0, -r, 0); // -r to maintain watertightness
    glRotated(90, 0, 0, 1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, metallic);
    glBegin(GL_QUAD_STRIP);
    for (i = 0; i < n_lengths; i++) {
        for (th = 360; th >= 0; th -= d_th) {
            glNormal3f(Cos(th), Sin(th), Sin(th));
            glTexCoord2d(Sin(th), 1);  glVertex3d(r * Cos(th), i * length_step + length_step, r * Sin(th));
            glTexCoord2d(Sin(th), 0);  glVertex3d(r * Cos(th), i * length_step              , r * Sin(th));
        }
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);

    glTranslated(0, n_lengths * length_step, 0);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, metallic);
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 1, 0);
    glTexCoord2d(0.5, 0.5);
    glVertex3f(0, 0.1, 0);
    for (th = 360; th >= 0; th -= d_th) {
        // add a bit of outwardness to the edge lighting so it looks a little less flat
        glNormal3f(Cos(th), 1, Sin(th));
        glTexCoord2d(Cos(th), Sin(th));  glVertex3f(r * Cos(th), 0, r * Sin(th));
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    /********** Support bars **********/
    r /= 2; // quarter radius for the support bars
    n_lengths = 7;
    glPushMatrix();
    glTranslated(0, 6 * length_step + 0.2, 0);
    glRotated(60, 0, 0, 1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, metallic);
    glBegin(GL_QUAD_STRIP);
    for (i = 0; i < n_lengths; i++) {
        for (th = 360; th >= 0; th -= d_th) {
            glNormal3f(Cos(th), Sin(th), Sin(th));
            glTexCoord2d(Sin(th), 1);  glVertex3d(r * Cos(th), i * length_step + length_step, r * Sin(th));
            glTexCoord2d(Sin(th), 0);  glVertex3d(r * Cos(th), i * length_step              , r * Sin(th));
        }
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    r /= 2; // eight radius for the vertical support bars
    n_lengths = 2;
    glPushMatrix();
    glTranslated(-2.3 * length_step, 7.4 * length_step + 0.2, 0);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, metallic);
    glBegin(GL_QUAD_STRIP);
    for (i = 0; i < n_lengths; i++) {
        for (th = 360; th >= 0; th -= d_th) {
            glNormal3f(Cos(th), Sin(th), Sin(th));
            glTexCoord2d(Sin(th), 1);  glVertex3d(r * Cos(th), i * length_step + length_step, r * Sin(th));
            glTexCoord2d(Sin(th), 0);  glVertex3d(r * Cos(th), i * length_step              , r * Sin(th));
        }
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);

    length_step = 0.4;
    glTranslated(-0.5, 0.21, 0);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, metallic);
    glBegin(GL_QUAD_STRIP);
    for (i = 0; i < n_lengths; i++) {
        for (th = 360; th >= 0; th -= d_th) {
            glNormal3f(Cos(th), Sin(th), Sin(th));
            glTexCoord2d(Sin(th), 1);  glVertex3d(r * Cos(th), i * length_step + length_step, r * Sin(th));
            glTexCoord2d(Sin(th), 0);  glVertex3d(r * Cos(th), i * length_step              , r * Sin(th));
        }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    length_step = 0.5;
    n_lengths = 8;
    stopLight(-(n_lengths + 6) * length_step, n_lengths * length_step, -0.13, 90, 1);

}

/*
 * Lays n_stretches stretches of 2m square sidewalk
 */
static void laySidewalk(int n_stretches) {
    const double sidewalk_width = 2.0; // 2m ~= 6 feet
    const double sidewalk_height = 0.2; // 0.2m is apparently a standard sidewalk height
    const double edge_width = 0.1;

    glColor3fv(concrete);
    set_material(ASPHALT);
    glPushMatrix();
    glTranslated(0, sidewalk_height, edge_width);
    glRotated(90, 0, 0, 1);
    drawRectangle(sidewalk_width - 2 * edge_width, sidewalk_width * n_stretches,
                  3, 3,
                  1, n_stretches,
                  pavement, 1);
    glPopMatrix();

    // draw soft corners as a 45-degree-rotated edge before the fully vertical edge
    // right soft corner
    glPushMatrix();
    glTranslated(0, sidewalk_height / 2, 0);
    glRotated(90, 0, 0, 1); // draw_rectangle is vertical by default, so make horizontal
    glRotated(45, 0, 1, 0); // the "soft" part of the corner

    // sqrt(2) * edge_width to bring the corner all the way to road w/45 degree angle
    drawRectangle(sqrt(2) * edge_width, sidewalk_width * n_stretches,
                  3, 3,
                  0.2, n_stretches,
                  pavement, 1);
    glRotated(45, 0, 1, 0); // additional rotation to vertical
    glTranslated(0, 0,-sqrt(2) * edge_width * Sin(45));
    drawRectangle(edge_width, sidewalk_width * n_stretches,
                  1, 1,
                  0.2, n_stretches,
                  pavement, 1);
    glPopMatrix();

    // left soft corner
    glPushMatrix();
    glTranslated(0, sidewalk_height, sidewalk_width - edge_width);
    glRotated(90, 0, 0, 1); // draw_rectangle is vertical by default, so make horizontal
    glRotated(-45, 0, 1, 0); // the "soft" part of the corner
    // sqrt(2) * edge_width to bring the corner all the way to road w/45 degree angle
    drawRectangle(sqrt(2) * edge_width, sidewalk_width * n_stretches,
                  1, 1,
                  0.2, n_stretches,
                  pavement, 1);
    glRotated(-45, 0, 1, 0); // additional rotation to vertical
    glTranslated(edge_width, 0, sqrt(2) * edge_width * Sin(45));
    drawRectangle(edge_width, sidewalk_width * n_stretches,
                  1, 1,
                  0.2, n_stretches,
                  pavement, 1);
    glPopMatrix();
}

/*
 * Lays road and sidewalk (on either side) in 12-meter stretches
 */
static void layRoad(int n_stretches, unsigned int texture) {
    const double sidewalk_width = 2.0; // 2m ~= 6 feet
    const double lane_width = 3.7; // might not even need this
    const double road_width = lane_width * 2;
    const double tex_length = 12; // relative length of road texture

    glColor3fv(concrete); // a good color for both road and sidewalk

    set_material(ASPHALT);
    glPushMatrix();
    glRotated(90, 0, 0, 1);
    
    drawRectangle(road_width, tex_length * n_stretches,
                  5, 10 * n_stretches,
                  1, n_stretches,
                  texture, 1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0, 0, -sidewalk_width);
    // draw as many sidewalk blocks as fill up the full road length
    laySidewalk((tex_length / sidewalk_width) * n_stretches);
    glTranslated(0, 0, sidewalk_width + road_width);
    laySidewalk((tex_length / sidewalk_width) * n_stretches);
    glPopMatrix();
}

/*
 * Draws a pre-fab fence
 */
static void drawFence(int n_stretches, float scale) {
    const float fence_length = 95;
    int i;

    set_material(WOOD);

    glPushMatrix();
    glTranslatef(0, 1.5, 0); // final model is apparently approx. 3 meters tall
    glScalef(scale,scale,scale);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, fence);

    for (i = 0; i < n_stretches; i++) {
        glTranslated(fence_length, 0, 0);
        glCallList(fence_obj);
    }

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

/*
 * Draws a n_quads vertical patches of grass across an area with size specified
 * by width and length.
 */
static void drawGrass(float width, float length, int n_quads, float scale) {
    // place quads at an even rate along both axes
    float norm = width + length;
    int n_quads_wid = (width / norm) * n_quads;
    int n_quads_len = (length / norm) * n_quads;
    float wid_step = width / (float) n_quads_wid;
    float len_step = length / (float) n_quads_len;

    int i, j;
    double corner = sqrt(2) / 2;
    // create one display list
    GLuint index = glGenLists(1);

    glNewList(index, GL_COMPILE);
        glBegin(GL_QUADS);
        glNormal3f(-1, 0, 1);
        glTexCoord2f(0, 0); glVertex3f(-corner, 0,-corner);
        glTexCoord2f(1, 0); glVertex3f( corner, 0, corner);
        glTexCoord2f(1, 1); glVertex3f( corner, 1, corner);
        glTexCoord2f(0, 1); glVertex3f(-corner, 1,-corner);

        glNormal3f(1, 0, -1);
        glTexCoord2f(0, 0); glVertex3f(-corner, 0,-corner);
        glTexCoord2f(0, 1); glVertex3f(-corner, 1,-corner);
        glTexCoord2f(1, 1); glVertex3f( corner, 1, corner);
        glTexCoord2f(1, 0); glVertex3f( corner, 0, corner);

        glNormal3f(1, 0, 1);
        glTexCoord2f(1, 0); glVertex3f(-corner, 0, corner);
        glTexCoord2f(0, 0); glVertex3f( corner, 0,-corner);
        glTexCoord2f(0, 1); glVertex3f( corner, 1,-corner);
        glTexCoord2f(1, 1); glVertex3f(-corner, 1, corner);

        glNormal3f(-1, 0, -1);
        glTexCoord2f(1, 0); glVertex3f(-corner, 0, corner);
        glTexCoord2f(1, 1); glVertex3f(-corner, 1, corner);
        glTexCoord2f(0, 1); glVertex3f( corner, 1,-corner);
        glTexCoord2f(0, 0); glVertex3f( corner, 0,-corner);
        glEnd();
    glEndList();

    set_material(GRASS);
    glUseProgram(grass_shader);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, grass);
    glUniform1i(glGetUniformLocation(grass_shader, "grass"), 0);

    for (i = 0; i < n_quads_wid; i++) {
        for (j = 0; j < n_quads_len; j++) {
            glPushMatrix();
            glTranslated(-width + i * wid_step, 0, -length + j * len_step);
            // rotate at random to make it look more convincing
            glRotated(grass_rot[i * n_quads_len + j], 0, 1, 0); // just like flattening a matrix
            glScalef(scale, scale, scale);
            glCallList(index);
            glPopMatrix();
        }
    }

    glDisable(GL_TEXTURE_2D);
    glDeleteLists(index, 1);
    glUseProgram(0);
}

/*
 *  Draw paper airplane
 *    at (x,y,z)
 *    nose towards (dx,dy,dz)
 *    up towards (ux,uy,uz)
 */
static void paperPlane()
{
    int i;
    // Dimensions used to size airplane
    const double wid=0.05;
    const double nose=+0.50;
    const double wing=+0.50;
    const double tail=-0.50;
    const double thic= 0.001;

    glColor3f(0.90, 0.90, 1.0);
    set_material(PLASTER); // not plaster, obviously, but it looks nice anyway
    /*
    * Below, pos_neg and extensive use of i and !i is to ensure that the wing
    * is always drawn counterclockwise
    */
    for (i = 0; i < 2; i++) {
        int pos_neg = 1 - i * 2; // 1 on the first iteration, -1 on the second
        int neg_pos =  -pos_neg; // -1 first iter, 1 second
        /********** "Fuselage" (v-shaped prism) **********/
        glBegin(GL_QUADS);
        // Topside
        glNormal3f(0, 1, neg_pos);
        glVertex3d(i * nose + !i * tail, wid, pos_neg * wid);
        glVertex3d(i * tail + !i * nose, wid, pos_neg * wid);
        glVertex3d(i * tail + !i * nose, 0  , 0            );
        glVertex3d(i * nose + !i * tail, 0  , 0            );
        // Underside
        glNormal3f(0,-1, pos_neg);
        glVertex3d(i * nose + !i * tail, wid - thic, pos_neg * wid);
        glVertex3d(i * nose + !i * tail,     - thic, 0            );
        glVertex3d(i * tail + !i * nose,     - thic, 0            );
        glVertex3d(i * tail + !i * nose, wid - thic, pos_neg * wid);
        glEnd();
        /********** Wings **********/
        glBegin(GL_TRIANGLES);
        // topside
        glNormal3f(0, 1, 0);
        glVertex3d(nose, wid, pos_neg * wid);
        glVertex3d(tail, wid, pos_neg * (!i * wid +  i * wing));
        glVertex3d(tail, wid, pos_neg * ( i * wid + !i * wing));

        // underside
        glNormal3f(0,-1, 0);
        glVertex3d(nose, wid - thic, pos_neg * wid);
        glVertex3d(tail, wid - thic, pos_neg * ( i * wid + !i * wing));
        glVertex3d(tail, wid - thic, pos_neg * (!i * wid +  i * wing));
        glEnd();
    }

    /********** Making everything watertight **********/
    glBegin(GL_QUADS);
    for (i = 0; i < 2; i++) {
        int pos_neg = 1 - i * 2; // 1 on the first iteration, -1 on the second
        // fuselage front face
        glNormal3f(1,0,0);
        glVertex3d(nose, wid             , pos_neg * wid);
        glVertex3d(nose,!i * (wid - thic),!i * wid      );
        glVertex3d(nose, - thic          , 0            );
        glVertex3d(nose, i * (wid - thic), i * (-wid)   );
        // fuselage back face
        glNormal3f(-1,0,0);
        glVertex3d(tail, wid             , pos_neg * wid);
        glVertex3d(tail, i * (wid - thic), i * (-wid)   );
        glVertex3d(tail, - thic          , 0            );
        glVertex3d(tail,!i * (wid - thic),!i * wid      );

        // wing back face
        glNormal3f(-1,0, 0);
        glVertex3d(tail, wid       , pos_neg * ( i * wing + !i * wid));
        glVertex3d(tail, wid - thic, pos_neg * ( i * wing + !i * wid));
        glVertex3d(tail, wid - thic, pos_neg * (!i * wing +  i * wid));
        glVertex3d(tail, wid       , pos_neg * (!i * wing +  i * wid));
        // wing diagonal face

        glNormal3f(i + !i / (wid - wing), 0, !i + i / (wid - wing));
        glVertex3d(tail, wid -  i * thic, pos_neg * wing);
        glVertex3d(tail, wid - !i * thic, pos_neg * wing);
        glVertex3d(nose, wid - !i * thic, pos_neg * wid );
        glVertex3d(nose, wid -  i * thic, pos_neg * wid );
    }
    glEnd();
}

/*
 * Convenience function for skySphere, transforms th, ph into texture coords
 */
static void skyTexCoord(int th, int ph) {
    ph = (ph + 90) / 2; // convert range to 0-90
    
    glTexCoord2f(0.5 + 0.5 * Cos(ph) * Cos(th),
                 0.5 + 0.5 * Cos(ph) * Sin(th));
}

/*
 * Draw a sphere surrounding the scene to apply a sky texture to
 */
static void skySphere(double r, unsigned int texture) {
    int th;
    int ph;
    int dTh = 5;
    int dPh = 5;

    if (zh < 180) {
        set_material(DAY);
    } else {
        set_material(NIGHT);
    }

    glColor3f(0.5 + 0.4 * Sin(zh), 0.5 + 0.4 * Sin(zh), 0.5 + 0.5 * Sin(zh)); // doesn't effect grass I don't think

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);

    // top
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, -1, 0);
    skyTexCoord(0, 90);
    glVertex3f(0, r, 0);
    for (th = 0; th <= 360; th += dTh) {
        ph = 90 - dPh;
        glNormal3f(-r * Cos(th) * Cos(ph),
                   -r * Sin(ph),
                   -r * Sin(th) * Cos(ph));
        skyTexCoord(th, ph);
        glVertex3f(r * Cos(th) * Cos(ph),
                   r * Sin(ph),
                   r * Sin(th) * Cos(ph));
    }
    glEnd();

    // main sphere
    glBegin(GL_QUAD_STRIP);
    for (ph = -90; ph < (90 - dPh); ph += dPh) {
        for (th = 360; th >= 0; th -= dTh) {
            glNormal3f(-r * Cos(th) * Cos(ph),
                       -r * Sin(ph),
                       -r * Sin(th) * Cos(ph));
            skyTexCoord(th, ph);
            glVertex3f(r * Cos(th) * Cos(ph),
                       r * Sin(ph),
                       r * Sin(th) * Cos(ph));

            glNormal3f(-r * Cos(th) * Cos(ph + dPh),
                       -r * Sin(ph + dPh),
                       -r * Sin(th) * Cos(ph + dPh));
            skyTexCoord(th, ph + dPh);
            glVertex3f(r * Cos(th) * Cos(ph + dPh),
                       r * Sin(ph + dPh),
                       r * Sin(th) * Cos(ph + dPh));
        }
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

/* 
 * Uses a mock-up of glider physics to update the position of the paper airplane
 * Definitions
 *  lift = lift_coeff * velocity^2 / 2
 *      coeff in this context also includes things like wing area, air density
 *  drag = coeff * velocity^2 / 2
 *      coeff in this context also includes reference area
 */
void physics_update() {
    double eff_t_diff = t_diff * time_scale;
    /* 
     * Find effective pitch by way of projecting forward airplane vector onto
     * the plane which defines the ground, then find angle between the forward
     * vector and its "ground component."
     */
    double projection[3];
    // stick the projection of the airplane forward vector & the ground into projection
    projectPlane3d(airplane[0], j_vec, projection);
    eff_pitch = findAngle3d(airplane[0], projection);
    eff_roll = findAngle3d(airplane[1], j_vec);

    double eff_lift_coeff = lift_coeff * fabs(Sin(eff_pitch * 2));
    lift = eff_lift_coeff * pow(velocity[0], 2) / 2;
    drag = drag_coeff * pow(velocity[0], 2) / 2;

    /*  
     * "Glide" is a mock force which nicely emulates gliding behavior
     *  with relatively little computation.
     * Honestly it does most of the heavy lifting here.
    */
    double glide = gravity * (0.1 + 0.9 * fabs(Sin(eff_pitch)));
    if (eff_pitch > 0 && airplane[0][1] > 0) {
        glide = -glide;
    }
    glide += push * push_val; push = 0;

    // acceleration due to drag is in the negative forward direction
    acceleration[0] = (glide -  drag - weight * airplane[0][1]) * eff_t_diff;
    acceleration[1] = (         lift - weight * airplane[1][1]) * eff_t_diff;
    acceleration[2] = ( 0            - weight * airplane[2][1]) * eff_t_diff;

    /* below,
     * velocity is the velocity of the airpane w.r.t. its own axis, so
     *  velocity[0] is the velocity of the airplane in the forward direction
     *  velocity[1] "  "   "        "  "   "        "  "   upward  "
     * and so on.
     */
    velocity[0] += acceleration[0];
    velocity[1] += acceleration[1];
    velocity[2] += acceleration[2];

    // keep velocity restricted to a sort of terminal velocity
    velocity[0] = limit(velocity[0], 4 * terminal_velocity);
    velocity[1] = limit(velocity[1], terminal_velocity);
    velocity[2] = limit(velocity[2], terminal_velocity);

    if (pos[1] > 0) { // allow crashing into ground
        pos[0] += (velocity[0] * airplane[0][0] + 
                   velocity[1] * airplane[1][0] +
                   velocity[2] * airplane[2][0]) * eff_t_diff;
        pos[1] += (velocity[0] * airplane[0][1] + 
                   velocity[1] * airplane[1][1] +
                   velocity[2] * airplane[2][1]) * eff_t_diff;
        pos[2] += (velocity[0] * airplane[0][2] + 
                   velocity[1] * airplane[1][2] +
                   velocity[2] * airplane[2][2]) * eff_t_diff;
    }
}

/*
 * Draws the objects which make up the scene (but not the ground or airplane)
 */
void scene() {
    /********** Skybox **********/
    glPushMatrix();
    glTranslated(0, -4 * dim, 0);
    skySphere(8 * dim, sky);
    glPopMatrix();

    /********** Fog **********/
    if (fog) {
        GLuint fogMode = GL_EXP2;
        GLfloat fogColor[4] = {0.55, 0.55, 0.7};
        glFogi(GL_FOG_MODE, fogMode);
        glFogfv(GL_FOG_COLOR, fogColor);
        glFogf(GL_FOG_DENSITY, 0.03);
        glHint(GL_FOG_HINT, GL_NICEST); // or GL_FASTEST
        glFogf(GL_FOG_START, 10.0);
        glFogf(GL_FOG_END, 4 * dim);
        glEnable(GL_FOG);
    } else {
        glDisable(GL_FOG);
    }


    /********** Ground and Scenery **********/
    /* below, most constants are experimentally determined */
    // draw fence along -z edge of map
    glPushMatrix();
    glTranslated(-dim - 2, 0, -dim);
    drawFence(15, 0.05);
    glPopMatrix();
    // draw fence along +x edge of map
    glPushMatrix();
    glTranslated(-dim, 0, -dim - 3);
    glRotated(-90, 0, 1, 0);
    drawFence(1, 0.05);
    glPopMatrix();
    // draw fence along -x edge of map
    glPushMatrix();
    glTranslated(dim, 0, -dim - 3);
    glRotated(-90, 0, 1, 0);
    drawFence(1, 0.05);
    glPopMatrix();
    // draw fence to enclose grass
    glPushMatrix();
    glTranslated(dim, 0, -22.5);
    glRotated(-90, 0, 1, 0);
    drawFence(8, 0.05);
    glPopMatrix();
    // opposite edge of map
    glPushMatrix();
    glTranslated(-dim, 0, -22.5);
    glRotated(-90, 0, 1, 0);
    drawFence(8, 0.05);
    glPopMatrix();
    // and finally the entire +z edge of map
    glPushMatrix();
    glTranslated(-dim - 2, 0, dim/2);
    drawFence(15, 0.05);
    glPopMatrix();

    glPushMatrix();
    glTranslated(dim, 0, -30);
    layRoad(3, road);
    glTranslated(-3 * 12, 0, 0); // translate by the 3 already-layed road lengths
    layRoad(1, xing);
    glTranslated(-12, 0, 0); // translate by 1 road length for crosswalk
    layRoad(2, road);
    glPopMatrix();

    // sidewalk leading up to apartment door
    glPushMatrix();
    glTranslated(-12, 0, -20.72);
    glRotated(90, 0, 1, 0);
    laySidewalk(10);
    glPopMatrix();
    // apartment door
    glPushMatrix();
    glTranslated(-11.85, 0.25, -2.1);
    glRotated(90, 0, 1, 0);
    rectangularPrism(3, 1.75, 0.05, 2, 3, 1, 1, door, 1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-6, 0, -21);
    glRotated(-90, 0, 1, 0);
    trafficLight();
    glPopMatrix();

    glPushMatrix();
    glTranslated(dim, 0, -2.1);
    drawGrass(dim * 1.25, 18, 64, 0.75);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-12, 0, -2.1);
    drawGrass(dim * 0.65, 18, 64, 0.75);
    glPopMatrix();

    glPushMatrix();
    glTranslated(dim, 0, -32);
    drawGrass(dim * 2, 3, 64, 2);
    glPopMatrix();

    /********** House **********/
    glPushMatrix();
    glTranslated(-11 * 2, 0, 3);
    drawApartment();
    glPopMatrix();
}

/*
 *  OpenGL (GLUT) calls this routine to display the scene
 */
void display()
{
    /********** Settings **********/
    //  Erase the window and the depth buffer
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    //  Enable Z-buffering in OpenGL
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    //  Undo previous transformations
    glLoadIdentity();
    glShadeModel(smooth ? GL_SMOOTH : GL_FLAT);

    float look_at[3]; // different from airplane pos array

    if (look_light) {
        look_at[0] = -6;
        look_at[1] = 4.5;
        look_at[2] = -27.75;
    } else if (look_grass) {
        look_at[0] = 0;
        look_at[1] = 1;
        look_at[2] = -dim + 2;
    } else if (look_window) {
        look_at[0] = -13.5;
        look_at[1] = 12.5;
        look_at[2] = -1;
    }
    
    if (look_light || look_grass || look_window) {
        double Ex = -0.1*dim*Sin(th)*Cos(ph);
        double Ey = +0.1*dim        *Sin(ph);
        double Ez = +0.1*dim*Cos(th)*Cos(ph);

        gluLookAt(look_at[0]-Ex,look_at[1]-Ey,look_at[2]-Ez , look_at[0],look_at[1],look_at[2] , 0,Cos(pitch),0);
    } else {
        // look at the airplane
        normalize3d(airplane[0]);
        // set ModelView matrix here so that lighting placement calculations are correct
        gluLookAt(pos[0] - airplane[0][0] * 5 , pos[1] - airplane[0][1] * 5 + 0.4 , pos[2] - airplane[0][2] * 5 ,
                  pos[0] , pos[1] , pos[2],
                  airplane[1][0] , airplane[1][1] , airplane[1][2]);
    }
    
    /********** Lighting **********/
    //  Light switch
    if (light)
    {
        //  Enable lighting with normalization
        glEnable(GL_LIGHTING);
        glEnable(GL_NORMALIZE);
        //  glColor sets ambient and diffuse color materials
        glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);

        /*
        * Light is turned off completely between 225 and 315 (when the sun
        * has dipped completely below the horizon), and its light is dampened
        * or partially dimmed on the intervals (180, 225) and (315, 360)
        */
        float dampening; // between 0 and 1, determines ambient light
        int input = 45;

        if (zh > 180 && zh < 225)
            input = 225 - zh; // 45 to 0
        else if (zh > 315 && zh < 360)
            input = 45 - (360 - zh); // 0 to 45
    
        dampening = Sin(2 * input);
        dampening = dampening < 0? 0 : dampening;

        float damped_one = 1.0 * dampening;
        float damped_point_one = 0.1 * dampening;

        //  Translate intensity to color vectors
        float Ambient[]   = {damped_point_one,damped_point_one,damped_point_one,damped_point_one};
        float Diffuse[]   = {damped_one,damped_one,damped_one,damped_one};
        float Specular[]  = {damped_point_one,damped_point_one,damped_point_one,damped_point_one};
        float white[]     = {1,1,1,1};

        float Position[] = {-3 * dim * Cos(zh),
                         3 * dim * Sin(zh),
                        -3 * dim * Sin(zh), 1};

        if (zh < 225 || zh > 315) {
            glEnable(GL_LIGHT0);
            glLightfv(GL_LIGHT0,GL_AMBIENT ,Ambient);
            glLightfv(GL_LIGHT0,GL_DIFFUSE ,Diffuse);
            glLightfv(GL_LIGHT0,GL_SPECULAR,Specular);
            glLightfv(GL_LIGHT0,GL_POSITION,Position);
            glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,32.0f);
            glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
            glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Diffuse);
        } else {
            glDisable(GL_LIGHT0);
        }

        // draw the sun
        ball(Position[0],Position[1],Position[2] , 8);
    }
    else
        glDisable(GL_LIGHTING);

    // discretize rotation, mouse doesn't really update nicely enough to store them
    if (roll != 0 || pitch != 0) {
        rotateVector3d(pitch, airplane[0], airplane[2]);
        crossProduct3d(airplane[0], airplane[2], airplane[1]);
        // perform cross product after each rotation to keep vectors orthogonal
        rotateVector3d(roll, airplane[1], airplane[0]);
        crossProduct3d(airplane[1], airplane[0], airplane[2]);

        pitch = 0;
        roll = 0;
    }
    

    /********** Physics **********/
    /* 
    * Update t values as close to their actual timing as possible,
    * makes the physics particularly accurate, but also makes the
    * airplane shake a little, presumably because of the difference
    * between when it is moved and when it is drawn.
    */
    t_diff = t - t_prev;
    physics_update(); // updates pos
    // update physics right after the effects are created
    t_prev = t;
    

    /********** PLANE **********/
    glPushMatrix();

    double airplane_col_major[16];
    glTranslated(pos[0], pos[1], pos[2]);

    // rotate so that plane displays in proper direction
    rotateVector3d(180, airplane[2], airplane[1]);
    
    // transform airplane to col major as openGL expects
    flattenToColMajor4d(airplane, airplane_col_major);
    glMultMatrixd(airplane_col_major);

    // rotate so that plane displays in proper direction
    rotateVector3d(180, airplane[2], airplane[1]);
    
    glScaled(1, 1, 1);
    paperPlane();

    glPopMatrix();

    /********** Draw the entire scene **********/
    //  Draw floor
    set_material(GRASS);
    glBindTexture(GL_TEXTURE_2D, sod);
    glEnable(GL_TEXTURE_2D);
    glNormal3f(0,1,0);
    int i, j;
    for (j=-dim;j<dim/2;j++)
    {
        glBegin(GL_QUAD_STRIP);
        for (i=-dim;i<=dim;i++)
        {
            glTexCoord2f(i,j);   glVertex3f(i,yFloor,j);
            glTexCoord2f(i,j+1); glVertex3f(i,yFloor,j+1);
        }
        glEnd();
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_TEXTURE_2D);

    scene();

    if (print) {
        glColor3d(1,1,1);
        //  Display parameters
        glWindowPos2i(5,45);
        Print("velocity[0]=%.6f   velocity[1]=%.6f   velocity[2]=%.6f", velocity[0], velocity[1], velocity[2]);
        glWindowPos2i(5,25);
        Print("zh=%i   time_scale=%.6f", zh, time_scale);
        glWindowPos2i(5,5);
        Print("eff_pitch=%i   eff_roll=%i\n", eff_pitch, eff_roll);
    }
    //  Render the scene and make it visible
    ErrCheck("display");
    glFlush();
    glutSwapBuffers();
}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
void special(int key,int x,int y)
{
    //  Right arrow key - increase angle by 5 degrees
    if (key == GLUT_KEY_RIGHT) {
        th += 5;
        roll += 5;
    }
    //  Left arrow key - decrease angle by 5 degrees
    else if (key == GLUT_KEY_LEFT) {
        th -= 5;
        roll -= 5;
    }
    //  Up arrow key - increase elevation by 5 degrees
    else if (key == GLUT_KEY_UP) {
        ph += 5;
        pitch += 5;
    }
    //  Down arrow key - decrease elevation by 5 degrees
    else if (key == GLUT_KEY_DOWN) {
        ph -= 5;
        pitch -= 5;
    }
    //  Keep angles to +/-360 degrees
    th %= 360;
    ph %= 360;
    roll %= 360;
    pitch %= 360;
    //  Update projection
    Project(fov,asp,dim);
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
void key(unsigned char ch,int x,int y)
{
    //  Exit on ESC
    if (ch == 27)
        exit(0);
    /********** Light movement **********/
    //  Toggle light
    else if (ch == 'l')
        light = 1-light;
    else if (ch == 'L')
        smooth = !smooth;
    else if (ch == 'f' || ch == 'F')
        fog = !fog;
    else if (ch == 'q' || ch == 'Q')
        time_scale /= 2;
    else if (ch == 'e' || ch == 'E')
        time_scale *= 2;
    else if (ch == 's' || ch == 'S')
        short_day = !short_day;
    else if (ch == 32) // space gives the plane a "push"
        push = !push;
    else if (ch == '0')
        look_light = look_grass = look_window = 0;
    else if (ch == '1') {
        th = 90;
        look_grass = look_window = 0;
        look_light = 1;
    } else if (ch == '2') {
        look_light = look_window = 0;
        look_grass = 1;
    } else if (ch == '3') {
        th = 180;
        ph = -15;
        look_light = look_grass = 0;
        look_window = 1;
    } else if (ch == 'p' || ch == 'P') {
        print = !print;
    } else if (ch == 'r' || ch == 'R') {
        // restart
        pos[0] = 0;
        pos[1] = 15;
        pos[2] = -15;
    }

    //  Reproject
    Project(fov,asp,dim);
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

void mouse_look(int x, int y) {
    roll -= X - x;
    pitch += Y - y;

    roll %= 180;
    pitch %= 180;

    X = x;
    Y = y;

    // Try to do some kind of wrap so that mouse doesn't get stuck at edges
    int windowWidth = glutGet(GLUT_WINDOW_WIDTH) - 5;
    int windowHeight = glutGet(GLUT_WINDOW_HEIGHT) - 5;

    if (X < 0 || X > windowWidth || Y < 0 || Y > windowHeight) {
        glutWarpPointer(windowWidth / 2, windowHeight / 2);
        X = windowWidth / 2;
        Y = windowHeight / 2;
    }

    glutPostRedisplay();
}

/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width,int height)
{
    //  Ratio of the width to the height of the window
    asp = (height>0) ? (double)width/height : 1;
    //  Set the viewport to the entire window
    glViewport(0,0, width,height);
    //  Set projection
    Project(fov,asp,dim);
}

/*
 *  Read text file
 */
char* ReadText(char *file)
{
   int   n;
   char* buffer;
   //  Open file
   FILE* f = fopen(file,"rt");
   if (!f) Fatal("Cannot open text file %s\n",file);
   //  Seek to end to determine size, then rewind
   fseek(f,0,SEEK_END);
   n = ftell(f);
   rewind(f);
   //  Allocate memory for the whole file
   buffer = (char*)malloc(n+1);
   if (!buffer) Fatal("Cannot allocate %d bytes for text file %s\n",n+1,file);
   //  Snarf the file
   if (fread(buffer,n,1,f)!=1) Fatal("Cannot read %d bytes for text file %s\n",n,file);
   buffer[n] = 0;
   //  Close and return
   fclose(f);
   return buffer;
}

/*
 *  Print Shader Log
 */
void PrintShaderLog(int obj,char* file)
{
   int len=0;
   glGetShaderiv(obj,GL_INFO_LOG_LENGTH,&len);
   if (len>1)
   {
      int n=0;
      char* buffer = (char *)malloc(len);
      if (!buffer) Fatal("Cannot allocate %d bytes of text for shader log\n",len);
      glGetShaderInfoLog(obj,len,&n,buffer);
      fprintf(stderr,"%s:\n%s\n",file,buffer);
      free(buffer);
   }
   glGetShaderiv(obj,GL_COMPILE_STATUS,&len);
   if (!len) Fatal("Error compiling %s\n",file);
}

/*
 *  Print Program Log
 */
void PrintProgramLog(int obj)
{
   int len=0;
   glGetProgramiv(obj,GL_INFO_LOG_LENGTH,&len);
   if (len>1)
   {
      int n=0;
      char* buffer = (char *)malloc(len);
      if (!buffer) Fatal("Cannot allocate %d bytes of text for program log\n",len);
      glGetProgramInfoLog(obj,len,&n,buffer);
      fprintf(stderr,"%s\n",buffer);
   }
   glGetProgramiv(obj,GL_LINK_STATUS,&len);
   if (!len) Fatal("Error linking program\n");
}

/*
 *  Create Shader
 */
int CreateShader(GLenum type,char* file)
{
   //  Create the shader
   int shader = glCreateShader(type);
   //  Load source code from file
   char* source = ReadText(file);
   glShaderSource(shader,1,(const char**)&source,NULL);
   free(source);
   //  Compile the shader
   fprintf(stderr,"Compile %s\n",file);
   glCompileShader(shader);
   //  Check for errors
   PrintShaderLog(shader,file);
   //  Return name
   return shader;
}

/*
 *  Create Shader Program
 */
int CreateShaderProg(char* VertFile,char* FragFile)
{
   //  Create program
   int prog = glCreateProgram();
   //  Create and compile vertex shader
   int vert = CreateShader(GL_VERTEX_SHADER  ,VertFile);
   //  Create and compile fragment shader
   int frag = CreateShader(GL_FRAGMENT_SHADER,FragFile);
   //  Attach vertex shader
   glAttachShader(prog,vert);
   //  Attach fragment shader
   glAttachShader(prog,frag);
   //  Link program
   glLinkProgram(prog);
   //  Check for errors
   PrintProgramLog(prog);
   //  Return name
   return prog;
}

/*
 *  GLUT calls this routine when the window is resized
 */
void idle()
{
    //  Elapsed time in seconds
    t = glutGet(GLUT_ELAPSED_TIME)/1000.0;

    t_diff = t - t_prev;
    physics_update(); // updates pos
    // update physics right after the effects are created
    t_prev = t;
        
    // day lasts either 18 or 360 seconds (5.5 minutes)
    if (short_day)
        zh = fmod(20 * t,360.0);
    else
        zh = fmod(t,360.0);
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  Start up GLUT and tell it what to do
 */
int main(int argc,char* argv[])
{
    //  Initialize GLUT
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(1440,900);
    glutCreateWindow("Luke Schaack Simple Flight");

    //  Set callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutSpecialFunc(special);
    glutKeyboardFunc(key);
    glutPassiveMotionFunc(mouse_look);
    glutIdleFunc(idle);

    //  Load textures
    metallic = LoadTexBMP("./textures/brushed_metal.bmp");
    hardwood = LoadTexBMP("./textures/hardwood2.bmp");
    tabletop = LoadTexBMP("./textures/tabletop.bmp");
    pavement = LoadTexBMP("./textures/pavement.bmp");
    plaster = LoadTexBMP("./textures/plaster4.bmp");
    shingle = LoadTexBMP("./textures/shingle.bmp");
    plastic = LoadTexBMP("./textures/plastic.bmp");
    siding = LoadTexBMP("./textures/siding_white.bmp");
    paint = LoadTexBMP("./textures/white_paint.bmp");
    grass = LoadTexBMP("./textures/grass2_lres.bmp");
    slats = LoadTexBMP("./textures/wood_slats.bmp");
    fence = LoadTexBMP("./textures/old_fence.bmp");
    soda = LoadTexBMP("./textures/soda2.bmp");
    road = LoadTexBMP("./textures/road.bmp");
    door = LoadTexBMP("./textures/door.bmp");
    xing = LoadTexBMP("./textures/xing.bmp");
    sky = LoadTexBMP("./textures/skydome.bmp");
    sod = LoadTexBMP("./textures/grass2.bmp");
    lid = LoadTexBMP("./textures/soda-top.bmp");

    //  Load objects
    fence_obj = LoadOBJ("./objects/old_fence.obj");

    //  Load shaders
    grass_shader = CreateShaderProg("./shaders/grass.vert","./shaders/grass.frag");
    tex_loc = glGetUniformLocation(grass_shader, "grass");

    // get some random numbers for grass rotation
    int i;
    for (i = 0; i < n_grass; i++) {
        grass_rot[i] = rand() % 90; // get a random number in the interval [0, 90)
    }

    //  Pass control to GLUT so it can interact with the user
    ErrCheck("init");
    glutMainLoop();
    return 0;
}
