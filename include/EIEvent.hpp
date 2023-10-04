#pragma once


#include "../external/Nucleus/include/Nucleus.hpp"


class EIEvent
{
protected:

    const long unsigned int& m_event_id;
    double m_Q;
    double m_Delta;
    Nucleus& m_nucleus;

public:

    virtual void sample();

    EIEvent() = delete;
    EIEvent(const long unsigned int& event_id, double Q, double Delta, Nucleus& nucleus);
    EIEvent(const EIEvent& other);
    EIEvent& operator=(const EIEvent&) = delete;
    EIEvent& operator=(EIEvent&&) = delete;
    virtual ~EIEvent() = default;
};


class EIEventCoherent : public EIEvent
{
    double x1, x2;
    double y1, y2;

public:

    void sample() override;
    double get_A() const;

    EIEventCoherent() = delete;
    EIEventCoherent(const EIEvent& super);
    EIEventCoherent(const EIEventCoherent&) = delete;
    EIEventCoherent& operator=(const EIEventCoherent&) = delete;
    EIEventCoherent& operator=(EIEventCoherent&&) = delete;
    ~EIEventCoherent() = default;
};


class EIEventIncoherent : public EIEvent
{
    double x1, x2;
    double y1, y2;
    double xb1, xb2;
    double yb1, yb2;

public:

    void sample() override;
    double get_A() const;

    EIEventIncoherent() = delete;
    EIEventIncoherent(const EIEvent& super);
    EIEventIncoherent(const EIEventIncoherent&) = delete;
    EIEventIncoherent& operator=(const EIEventIncoherent&) = delete;
    EIEventIncoherent& operator=(EIEventIncoherent&&) = delete;
    ~EIEventIncoherent() = default;
};