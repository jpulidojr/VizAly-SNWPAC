/* wavelet/wavelet.c
 * 
 * Copyright (C) 2004, 2009 Ivo Alxneit
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Modified by Jesus Pulido */

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_wavelet.h>

#define RETURN_IF_NULL(x) if (!x) { return ; }

gsl_wavelet *
gsl_wavelet_alloc (const gsl_wavelet_type * T, size_t k)
{
  int status;

  gsl_wavelet *w = (gsl_wavelet *) malloc (sizeof (gsl_wavelet));

  if (w == NULL)
    {
	  printf("failed to allocate space for wavelet struct"); abort();
    };

  w->type = T;

  status = (T->init) (&(w->h1), &(w->g1), &(w->h2), &(w->g2),
                      &(w->nc), &(w->offset), k);

  if (status)
    {
      free (w);
	  printf("invalid wavelet member"); abort();
    }

  return w;
}

void
gsl_wavelet_free (gsl_wavelet * w)
{
  RETURN_IF_NULL (w);
  free (w);
}

const char *
gsl_wavelet_name (const gsl_wavelet * w)
{
  return w->type->name;
}


/* Let's not export this for now (BJG) */
#if 0
void
gsl_wavelet_print (const gsl_wavelet * w)
{
  size_t n = w->nc;
  size_t i;

  printf ("Wavelet type: %s\n", w->type->name);

  printf
    (" h1(%d):%12.8f   g1(%d):%12.8f       h2(%d):%12.8f   g2(%d):%12.8f\n",
     0, w->h1[0], 0, w->g1[0], 0, w->h2[0], 0, w->g2[0]);

  for (i = 1; i < (n < 10 ? n : 10); i++)
    {
      printf
        (" h1(%d):%12.8f   g1(%d):%12.8f       h2(%d):%12.8f   g2(%d):%12.8f\n",
         i, w->h1[i], i, w->g1[i], i, w->h2[i], i, w->g2[i]);
    }

  for (; i < n; i++)
    {
      printf
        ("h1(%d):%12.8f  g1(%d):%12.8f      h2(%d):%12.8f  g2(%d):%12.8f\n",
         i, w->h1[i], i, w->g1[i], i, w->h2[i], i, w->g2[i]);
    }
}
#endif

gsl_wavelet_workspace *
gsl_wavelet_workspace_alloc (size_t n)
{
  gsl_wavelet_workspace *work;

  if (n == 0)
    {
	  printf("length n must be positive integer"); abort();
    }

  work = (gsl_wavelet_workspace *) malloc (sizeof (gsl_wavelet_workspace));

  if (work == NULL)
    {
	  printf("failed to allocate struct"); abort();
    }

  work->n = n;
  work->scratch = (double *) malloc (n * sizeof (double));

  if (work->scratch == NULL)
    {
      /* error in constructor, prevent memory leak */
      free (work);
	  printf("failed to allocate scratch space"); abort();
    }

  return work;
}

void
gsl_wavelet_workspace_free (gsl_wavelet_workspace * work)
{
  RETURN_IF_NULL (work);
  /* release scratch space */
  free (work->scratch);
  work->scratch = NULL;
  free (work);
}
