document.addEventListener('DOMContentLoaded', function() {
    var bodyTag = document.getElementsByTagName('body')[0];

    var text = `
    <div class="navbar fixed-top navbar-color page-columns page-layout-article">
        <div id="quarto-margin-sidebar" class="sidebar margin-sidebar"></div>
        <div class="content">
            <div class="d-flex navbar-report align-items-end">
                <div class="me-auto" style="font-weight: 500">
                    <img style="height: 30px" src="./static/images/logo.png">
                </div>
                <div style="padding-left: 1rem;">
                    <div></div>
                </div>
            </div>
        </div>
    </div>`;

    bodyTag.insertAdjacentHTML("afterbegin", text);
});
