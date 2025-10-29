import pyncare as pm


def q_plot(ax, qfactor: pm.Qfactor):
    """Plots the q factor extraced and derived data, as a function of ψp."""
    psip_data = qfactor.psip_data()
    q_data = qfactor.q_data()
    q_data_derived = qfactor.q_data_derived()

    ax.scatter(psip_data, q_data, c="k", s=2, zorder=2, alpha=0.8, label="data points")
    ax.plot(psip_data, q_data, c="b", label=r"$q(\psi_p)$")
    ax.plot(psip_data, q_data_derived, c="r", label=r"$d\psi / d\psi_p$")

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$q(\psi_p)$")
    ax.set_title(r"$q(\psi_p)$")
    ax.grid(True)
    ax.margins(0)
    ax.legend()


def psi_plot(ax, qfactor: pm.Qfactor):
    """Plots the ψ extraced data, as a function of ψp."""
    psip_data = qfactor.psip_data()
    psi_data = qfactor.psi_data()

    ax.scatter(
        psip_data, psi_data, c="k", s=2, zorder=2, alpha=0.8, label="data points"
    )
    ax.plot(psip_data, psi_data, c="b", label=r"$\psi(\psi_p)$")

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$\psi(\psi_p)$")
    ax.set_title(r"$\psi(\psi_p)$")
    ax.grid(True)
    ax.margins(0)
    ax.legend()
