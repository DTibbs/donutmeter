/* Copyright (C) 2017 Rainmeter Project Developers
 *
 * This Source Code Form is subject to the terms of the GNU General Public
 * License; either version 2 of the License, or (at your option) any later
 * version. If a copy of the GPL was not distributed with this file, You can
 * obtain one at <https://www.gnu.org/licenses/gpl-2.0.html>. */

// Required for M_PI in cmath
#define _USE_MATH_DEFINES

#include <Windows.h>
#include <RainmeterAPI.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <cwctype>
#include <string>
#include <stdint.h>
#include <vector>

//
// See <Project root>/Skins/donutmeter/donutmeter.ini for Rainmeter Skin usage
//

struct DonutMeasure
{
    // How many chars to render per row
    uint32_t m_screen_width;
    // How many rows of chars to render
    uint32_t m_screen_height;
    // Cache pixel count (width * height) for calculating indexes later
    uint32_t m_pixel_count;
    // R1 - The girth of the donut
    double m_r1;
    // R2 - The girth + inner radius of donut
    double m_r2;
    // K1 - This will be calculated based on the screen width, with a ratio from edge of screen
    double m_k1;
    // K2 - Distance of donut's position from origin (eye is considered at origin)
    double m_k2;
    // phi_spacing - This will be calculated by dividing 2PI / number of segments
    //
    // Number of donut segments. phi_spacing = 2PI / segment count
    //        E.g. 8 segments == 2PI / 8 == PI / 4 == ~0.79 or 45 degrees
    //
    //                       Y-axis
    //                          |
    //                                           
    //                     N    V   -                  
    //                XK0Okd    V   R                   
    //             Xkolcccc;    V   1        //            
    //           Xklc::;;cc;    V - -      // 
    //          Ko:c:coc;;;,    V |      // 
    //         Nkc;:;:l:,,:;    V R    // 
    //        NOlcolcccc,';l    V 2  // 
    //        kol:c:::;,':ON      |          
    //       Xl;c;;;,:l;;0        |                   
    //      ------------------    -   ===============    ------- X-axis
    //       Xdoo::c;','cX             
    //       Xd::;;:;''';O             
    //        k::c:;;;,'':kN         \\ 
    //        Nkoccoc:;''',l    V      \\ 
    //         Nkl::,cl,,;;:    V        \\ 
    //          Xxcloc:;:c:,    V          \\ 
    //           NOdl:cccc;;    V            \\ 
    //        ^    NOdo:lxc;    V                      
    //         \      XOxdlc    V 
    //          \         NX    V 
    //           - etc.         
    //
    double m_phi_spacing;
    // thets_spacing - This will be calculated by dividing 2PI / number of segment points;
    //
    //  Number of donut segment points. theta_spacing = 2PI / segment point count
    //        E.g. 8 segments == 2PI / 8 == PI / 4 == ~0.79 or 45 degrees
    //
    //                    Y-axis
    //                       |                        
    //                       _                         
    //          NNNNNNNNN   |O|                  
    //              NNNNN    - _              
    //           NXN  NN     | |         _             
    //         NXN   NN      | |        |O|      
    //        NXN   N        | R         -             
    //       XX              | 1      
    //      NX               | |                       
    //      XN               | |                       
    //     NXN               | |                _     
    //---------------------------------------  |O|  ---- Z-axis
    //      XN               |                  -      
    //      XX               |                         
    //      NXN              |                        
    //       NXN             |                          
    //        NXN            |           _             
    //          XXN          |          |O|            
    //            NNNN       |           -             
    //     ^        NNNNNNN  _                         
    //      \           NNN |O|                        
    //       \               -                        
    //        - etc.         |                        
    //                       |                        
    //
    double m_theta_spacing;
    // The array of text to display rendered donut
    std::wstring m_output;
    // The depth buffer to calculate luminance of donut
    std::vector<double>  m_zbuffer;
    // Keep a timer to calculate rotations
    time_t m_timer;

    double m_A, m_B;

    // Cache value of 2 * M_PI for reuse
    static const double two_pi;

    // Cache luminance ascii characters.
    //  TODO: Make this configurable by skin?
    static const std::string lum_chars;

    DonutMeasure()
        : m_screen_width(0)
        , m_screen_height(0)
        , m_pixel_count(0)
        , m_r1(0)
        , m_r2(0)
        , m_k1(0)
        , m_k2(0)
        , m_phi_spacing(0)
        , m_theta_spacing(0) 
        , m_output()
        , m_zbuffer()
        , m_timer(std::time(nullptr))
        , m_A()
        , m_B() {}

    void render_frame(double A, double B);
};

// Define value of two_pi here
const double DonutMeasure::two_pi = 2.0 * M_PI;
// Define lum_chars here
const std::string DonutMeasure::lum_chars = ".,-~:;=!*#$@";

void DonutMeasure::render_frame(double A, double B)
{
    // precomputing sines/cosines
    double cosA = cos(A), sinA = sin(A);
    double cosB = cos(B), sinB = sin(B);

    // theta is for generating the points for each segment of the donut
    for(double theta = 0; theta < DonutMeasure::two_pi; theta += m_theta_spacing)
    {
        // precomppute sin/cos of theta
        double costheta = cos(theta), sintheta = sin(theta);
        // phi is for generating each segment of the donut
        for(double phi = 0; phi < DonutMeasure::two_pi; phi += m_phi_spacing)
        {
            // precompute sin/cos of phi
            double cosphi = cos(phi), sinphi = sin(phi);

            // x,y of circle center to generate points at
            double circlex = m_r2 + m_r1 * costheta;
            double circley = m_r1 * sintheta;

            // Calculate x,y,z with rotation
            double x = circlex * (cosB * cosphi + sinA * sinB * sinphi) - circley * cosA * sinB;
            double y = circlex * (sinB * cosphi - sinA * cosB * sinphi) - circley * cosA * cosB;
            double z = m_k2 + cosA * circlex * sinphi * circley * sinA;
            double invz = 1 / z;

            // Project to screen (centered), inverting y (2D y axis is opposite of 3D's y axis)
            int xp = (int)(m_screen_width / 2 + m_k1 * invz * x);
            int yp = (int)(m_screen_height / 2 - m_k1 * invz * y);

            // Calculate luminance (Light is at 0, 1, -1, which is baked into this formula.
            //  TODO: support changing the light's position?
            double lum = cosphi * costheta * sinB - cosA * costheta * sinphi - sinA * sintheta + cosB * (cosA * sintheta - costheta * sinA * sinphi);
            // Luminance ranges from -sqrt(2) to +sqrt(2). > 0 is facing the screen, so only use those that are facing that way.
            if(lum > 0)
            {
                int buf_idx = (yp * m_screen_width) + xp;
                
                // This shouldn't happen...
                if(buf_idx > m_pixel_count)
                {
                    continue;
                }

                //Test against zbuffer. Larger invz is closer to the screen.
                if(invz > m_zbuffer[buf_idx])
                {
                    m_zbuffer[buf_idx] = invz;
                    // Convert Luminance to a 0..11 scale by multiplying by 8 (e.g. 8*sqrt(2) == 11)
                    int lum_idx = lum * 8;
                    // This shouldn't happen...
                    if(lum_idx > 11)
                    {
                        continue;
                    }
                    m_output[buf_idx] = lum_chars[lum_idx];
                }
            }
        }
    }
}

PLUGIN_EXPORT void Initialize(void** data, void* rm)
{
    DonutMeasure* measure = new DonutMeasure;
    *data = measure;
}

PLUGIN_EXPORT void Reload(void* data, void* rm, double* maxValue)
{
    DonutMeasure* measure = (DonutMeasure*)data;

    // Load values from config on reload
    measure->m_screen_width = RmReadInt(rm, L"ScreenWidth", 100);
    measure->m_screen_height = RmReadInt(rm, L"ScreenHeight", 100);
    measure->m_pixel_count = measure->m_screen_width * measure->m_screen_height;
    measure->m_r1 = RmReadDouble(rm, L"DonutGirthRadius", 1.0);
    measure->m_r2 = measure->m_r1 + RmReadDouble(rm, L"DonutHoleRadius", 1.0);
    measure->m_k2 = RmReadDouble(rm, L"DonutZDistance", 5.0);
    measure->m_phi_spacing = DonutMeasure::two_pi / RmReadDouble(rm, L"DonutSegmentCount", 25);
    measure->m_theta_spacing = DonutMeasure::two_pi / RmReadDouble(rm, L"DonutSegmentPointCount", 12);
    //  Maximum x-distance of donut is m_r1 + m_r2 when z = 0. Using a 3/8ths displacement.
    measure->m_k1 = measure->m_screen_width * measure->m_k2 * 3 / (8 * (measure->m_r1 + measure->m_r2));

    measure->m_A = 0.0;
    measure->m_B = 0.0;

    // Clear output and zbuffer
    measure->m_output.assign(measure->m_pixel_count, ' ');
    measure->m_zbuffer.assign(measure->m_pixel_count, 0.0);
    // Put newlines at every 'edge'
    for(uint32_t row = 0; row < measure->m_screen_height; ++row)
    {
        measure->m_output[(measure->m_screen_width - 1) + (measure->m_screen_width * row)] = '\n';
    }
}

PLUGIN_EXPORT double Update(void* data)
{
    DonutMeasure* measure = (DonutMeasure*)data;

    // Clear output and zbuffer
    measure->m_output.assign(measure->m_pixel_count, ' ');
    measure->m_zbuffer.assign(measure->m_pixel_count, 0.0);
    // Put newlines at every 'edge'
    for(uint32_t row = 0; row < measure->m_screen_height; ++row)
    {
        measure->m_output[(measure->m_screen_width - 1) + (measure->m_screen_width * row)] = '\n';
    }

    measure->m_A += 0.04;
    measure->m_B += 0.02;
    measure->render_frame(measure->m_A, measure->m_B);
    //measure->m_output.assign(                        \
    //  L"// Obfuscated donut                       \n \
    //              k;double sin()                  \n \
    //           ,cos();main(){float A=             \n \
    //         0,B=0,i,j,z[1760];char b[            \n \
    //       1760];printf(\"\\x1b[2J\");for(;;      \n \
    //    ){memset(b,32,1760);memset(z,0,7040)      \n \
    //    ;for(j=0;6.28>j;j+=0.07)for(i=0;6.28      \n \
    //   >i;i+=0.02){float c=sin(i),d=cos(j),e=     \n \
    //   sin(A),f=sin(j),g=cos(A),h=d+2,D=1/(c*     \n \
    //   h*e+f*g+5),l=cos      (i),m=cos(B),n=s\\   \n \
    //  in(B),t=c*h*g-f*        e;int x=40+30*D*    \n \
    //  (l*h*m-t*n),y=            12+15*D*(l*h*n    \n \
    //  +t*m),o=x+80*y,          N=8*((f*e-c*d*g    \n \
    //   )*m-c*d*e-f*g-l        *d*n);if(22>y&&     \n \
    //   y>0&&x>0&&80>x&&D>z[o]){z[o]=D;;;b[o]=     \n \
    //   \".,-~:;=!*#$@\"[N>0?N:0];}}/*#****!!-*/   \n \
    //    printf(\"\\x1b[H\");for(k=0;1761>k;k++)   \n \
    //     putchar(k%80?b[k]:10);A+=0.04;B+=        \n \
    //       0.02;}}/*****####*******!!=;:~         \n \
    //         ~::==!!!**********!!!==::-           \n \
    //           .,~~;;;========;;;:~-.             \n \
    //               ..,--------,*/                 \n ");


    return 0.0;
}

PLUGIN_EXPORT LPCWSTR GetString(void* data)
{
    DonutMeasure* measure = (DonutMeasure*)data;
    return  measure->m_output.c_str();
}

PLUGIN_EXPORT void Finalize(void* data)
{
    DonutMeasure* measure = (DonutMeasure*)data;
    delete measure;
}
