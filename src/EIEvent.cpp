#include <random>

#include "../include/EIEvent.hpp"
#include "../include/utilities.hpp"
#include "../include/Coherent.hpp"
#include "../include/Incoherent.hpp"


void EIEvent::sample()
{
    m_nucleus.sample_pos();
}


EIEvent::EIEvent (const long unsigned int& event_id, double Q, double Delta, Nucleus& nucleus)
    : m_event_id(event_id), m_Q(Q), m_Delta(Delta), m_nucleus(nucleus)
{
}


EIEvent::EIEvent (const EIEvent& other)
    : m_event_id(other.m_event_id), m_Q(other.m_Q), m_Delta(other.m_Delta), m_nucleus(other.m_nucleus)
{
}


void EIEventCoherent::sample()
{
    x1 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
    x2 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());

    y1 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
    y2 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
}


double EIEventCoherent::get_A () const
{
    double A_sum = 0.0;

    for (uint i=0; i<m_nucleus.get_atomic_num(); i++)
    {
        const double* pos;
        
        pos = m_nucleus.get_nucleon_pos(i);

        A_sum += Coherent::A_integrand_function((x1+y1)/2.0-pos[0], (x2+y2)/2.0-pos[1], x1-y1, x2-y2, m_Q, m_Delta);
    }

    return A_sum;
}


EIEventCoherent::EIEventCoherent (const EIEvent& super)
    : EIEvent(super)
{
    sample();
}


void EIEventIncoherent::sample()
{
    x1 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
    x2 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());

    y1 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
    y2 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());

    xb1 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
    xb2 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());

    yb1 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
    yb2 = rng(-2.0*m_nucleus.get_bulk_radius(), 2.0*m_nucleus.get_bulk_radius());
}


double EIEventIncoherent::get_A() const
{
    double A_sum = 0.0;

    for (uint i=0; i<m_nucleus.get_atomic_num(); i++)
    {
        const double* pos;
        
        pos = m_nucleus.get_nucleon_pos(i);

        A_sum += Incoherent::A_integrand_function((x1+y1)/2.0-pos[0], (x2+y2)/2.0-pos[1], x1-y1, x2-y2, (xb1+yb1)/2.0-pos[0], (xb2+yb2)/2.0-pos[1], xb1-yb1, xb2-yb2, m_Q, m_Delta);
    }

    return A_sum;
}


EIEventIncoherent::EIEventIncoherent (const EIEvent& super)
    : EIEvent(super)
{
    sample();
}