#ifndef SIM_H
#define SIM_H

class sim {
private:
    static sim *instance;
public:
    inline static void setInstance(sim *s) { sim::instance = s; }
    inline static sim *getInstance() { return sim::instance; }
    virtual double time() const = 0;
    virtual void op() = 0;
};

#endif /* SIM_H */

