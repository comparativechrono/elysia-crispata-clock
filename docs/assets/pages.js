/* Hover and focus tooltips for chart marks on the article pages. */
(function () {
  "use strict";
  const tip = document.createElement("div");
  tip.className = "tip-pop"; tip.setAttribute("role", "tooltip");
  document.body.append(tip);
  function show(text, x, y) {
    tip.textContent = text; tip.classList.add("on");
    const r = tip.getBoundingClientRect();
    let left = x + 14, top = y + 14;
    if (left + r.width > innerWidth - 8) left = x - r.width - 14;
    if (top + r.height > innerHeight - 8) top = y - r.height - 14;
    tip.style.left = Math.max(6, left) + "px"; tip.style.top = Math.max(6, top) + "px";
  }
  const hide = () => tip.classList.remove("on");
  document.querySelectorAll("[data-tip]").forEach((n) => {
    const t = n.querySelector("title"); if (t) t.remove();   // avoid a second, native tooltip
    n.addEventListener("mousemove", (e) => show(n.dataset.tip, e.clientX, e.clientY));
    n.addEventListener("mouseleave", hide);
    n.addEventListener("focus", () => { const r = n.getBoundingClientRect(); show(n.dataset.tip, r.right, r.top); });
    n.addEventListener("blur", hide);
  });
})();
